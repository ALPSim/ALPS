import gc
import tempfile
import weakref

import numpy as np
# Initialize the package DLL directory before importing a native consumer.
from pyalps import hdf5, ngs
import parameter_probe as native


def fails(function, *arguments):
    try:
        function(*arguments)
    except (RuntimeError, TypeError, ValueError, IndexError, OverflowError) as error:
        return str(error)
    raise AssertionError(f"{function.__name__}{arguments} did not fail")


empty = native.empty_vectors()
for key, kind in {"integer": "i", "real": "f", "boolean": "b"}.items():
    assert empty[key].size == 0 and empty[key].dtype.kind == kind
assert isinstance(empty["text"], list) and empty["text"] == []

# Values a program defines read back with its type; string vectors as lists.
parameters = native.string_parameters()
assert parameters["value"] == ["Sz"]
parameters["value"] = ["Magnetization", "Susceptibility"]
assert native.string_vector(parameters) == ["Magnetization", "Susceptibility"]

# Python values are supplied as text when the params are made, so the
# params do not follow later changes to the Python object.
values = np.array([1.0, 2.0])
parameters = ngs.params({"value": values})
values *= 3
assert native.vector(parameters) == [1.0, 2.0]
assert native.threaded_vector(parameters) == [1.0, 2.0]
clone = native.clone(parameters)
native.replace(parameters)
assert native.vector(parameters) == [5.0, 6.0]
assert native.vector(clone) == [1.0, 2.0]
parameters["value"] = [1.5, 2.5, 3.5]
assert native.vector(parameters) == [1.5, 2.5, 3.5]
parameters = ngs.params({"value": np.array([2 ** 40, 2 ** 40 + 1], dtype=np.int64)})
assert native.vector(parameters) == [2 ** 40, 2 ** 40 + 1]
assert native.integer_vector(parameters) == [2 ** 40, 2 ** 40 + 1]

# Integers are exact and range-checked by the type the program defines.
parameters = ngs.params({"value": 2 ** 53 + 1})
assert native.wide_integer(parameters) == 2 ** 53 + 1
assert "value" in fails(native.integer, parameters)
assert native.wide_integer(ngs.params({"value": np.array(2 ** 53 + 1, dtype=np.int64)})) == 2 ** 53 + 1
# A non-integral number for an integer is an error naming the parameter.
assert "value" in fails(native.integer, ngs.params({"value": 2.5}))
assert native.integer(ngs.params({"value": 24.0})) == 24

# Python and NumPy Booleans are supplied as 1 and 0.
for value in (True, np.bool_(True)):
    assert native.flag(ngs.params({"value": value})) is True
    assert native.integer(ngs.params({"value": value})) == 1

# Mixed lists convert element by element with the program's type.
for value, real, integers in (
    ([True, 2, 3], [1., 2., 3.], [1, 2, 3]),
    ([1, "2.5", 3.5], [1., 2.5, 3.5], None),
    ((1, "2", 3), [1., 2., 3.], [1, 2, 3]),
    ([2 ** 53 + 1, 2], [float(2 ** 53 + 1), 2.], [2 ** 53 + 1, 2]),
):
    parameters = ngs.params({"value": value})
    assert native.vector(parameters) == real
    if integers is None:
        assert "value" in fails(native.integer_vector, parameters)
    else:
        assert native.integer_vector(parameters) == integers

for value, expected in (
    (["Energy", "Stiffness"], "Energy,Stiffness"),
    ([1, "two", 3.0], "1,two,3"),
    ((True, False), "1,0"),
    (["", "middle", ""], ",middle,"),
    ([], ""),
):
    assert native.text(ngs.params({"value": value})) == expected
assert native.string_vector(ngs.params({"value": ["Energy", "Stiffness"]})) == ["Energy", "Stiffness"]
# params holds no complex numbers, None, nested lists or elements with commas
for value in (2 + 3j, [1, 2 + 3j], None, [None], [[1, 2]], ["a,b"], np.ones((2, 3))):
    try:
        ngs.params({"value": value})
    except TypeError as error:
        assert "'value'" in str(error), error
    else:
        raise AssertionError(f"{value!r} was accepted")

with tempfile.TemporaryDirectory() as directory:
    filename = directory + "/shared-archive.h5"
    with hdf5.archive(filename, "w") as archive:
        archive["initial"] = 1
    with hdf5.archive(filename, "r") as archive:
        # Opening a writer upgrades the existing read-only archive context.
        # A second SDK runtime copy has a separate registry and cannot do this.
        native.append_archive(filename)
        assert archive["initial"] == 1
        assert archive["from_native"] == 42

    # A group pyalps wrote from a dict loads as text (no ALPSCore attributes);
    # a nested dict is not a parameter and is an error naming it.
    with hdf5.archive(directory + "/parameters.h5", "w") as archive:
        archive["parameters"] = {"value": np.array([2.0, 4.0]), "L": 4}
        archive["nested"] = {"metadata": {"label": "test"}}
        archive.set_context("/parameters")
        native.load(clone, archive)
        assert native.threaded_vector(clone) == [2.0, 4.0]
        assert archive.context == "/parameters"
        archive.set_context("/nested")
        assert "'metadata'" in fails(native.load, clone, archive)

    # Lists of Booleans or mixed scalars are numbered HDF5 groups. More than
    # ten elements catches lexical order.
    for i, values in enumerate((
        [True, False], [True, 2, 3],
        [True] + list(range(1, 12)),
        [True, np.int8(-7), np.int64(2 ** 53 + 1)],
    )):
        with hdf5.archive(directory + f"/native-list-{i}.h5", "w") as archive:
            archive["parameters"] = {"value": values}
            archive.set_context("/parameters")
            native_parameters = native.empty_vectors()
            native.load(native_parameters, archive)
            assert native.vector(native_parameters) == [float(value) for value in values]
            assert native.integer_vector(native_parameters) == [int(value) for value in values]
            assert archive.context == "/parameters"
            # ALPSCore's own format keeps the values as supplied text.
            archive.set_context("/resaved")
            native.save(native_parameters, archive)
            restored = native.empty_vectors()
            native.load(restored, archive)
            assert native.vector(restored) == [float(value) for value in values]
            assert native.integer_vector(restored) == [int(value) for value in values]
            assert archive.context == "/resaved"

# Native destruction on a thread that started without the GIL holds no
# Python value, so the params never keep one alive.
value = np.ones(3)
reference = weakref.ref(value)
parameters = ngs.params({"value": value})
del value
gc.collect()
assert reference() is None
native.destroy_on_worker(parameters)
print("native parameter contracts: ok")
