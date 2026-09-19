import copy
import gc
import tempfile
import weakref

import numpy as np
import parameter_probe as native
from pyalps import hdf5, ngs

empty = native.empty_vectors()
for key, kind in {"integer": "i", "real": "f", "complex": "c", "boolean": "b", "text": "U"}.items():
    assert empty[key].size == 0 and empty[key].dtype.kind == kind

values = np.array([1.0, 2.0])
parameters = ngs.params({"value": values})
assert native.vector(parameters) == [1.0, 2.0]
values *= 3
assert native.vector(parameters) == [3.0, 6.0]
parameters["value"][0] = 8
assert native.threaded_vector(parameters) == [8.0, 6.0]
clone = native.clone(parameters)
parameters["value"][1] = 9
assert native.vector(clone) == [8.0, 9.0]
native.replace(parameters)
assert native.vector(parameters) == [5.0, 6.0]
assert native.vector(clone) == [8.0, 9.0]
parameters["value"][0] = 10
assert native.vector(parameters) == [10.0, 6.0]
parameters["value"] = [1.5, 2.5]
parameters["value"].append(3.5)
assert native.vector(parameters) == [1.5, 2.5, 3.5]
parameters["value"] = np.array([2 ** 40, 2 ** 40 + 1], dtype=np.int64)
assert native.vector(parameters) == [2 ** 40, 2 ** 40 + 1]
parameters["value"] = 2 ** 53 + 1
assert native.wide_integer(parameters) == 2 ** 53 + 1
try:
    native.integer(parameters)
except (RuntimeError, TypeError, ValueError, IndexError, OverflowError):
    pass
else:
    raise AssertionError("out-of-range conversion to native int must fail")
parameters["value"] = np.array(2 ** 53 + 1, dtype=np.int64)
assert native.wide_integer(parameters) == 2 ** 53 + 1

# The legacy reader converted each list element to the requested C++ type.
# A homogeneous intermediate vector rejects valid mixed inputs or loses an
# imaginary component before a complex-valued consumer can read it.
for value, real, complex_values, integers in (
    ([True, 2, 3.5], [1., 2., 3.5], [1+0j, 2+0j, 3.5+0j], [1, 2, 3]),
    ([1, "2.5", 3.5], [1., 2.5, 3.5], [1+0j, 2.5+0j, 3.5+0j], [1, 2, 3]),
    ([True, 2+3j, "4"], [1., 2., 4.], [1+0j, 2+3j, 4+0j], [1, 2, 4]),
    ((1, "2", 3), [1., 2., 3.], [1+0j, 2+0j, 3+0j], [1, 2, 3]),
    ([2**53+1, 2], [float(2**53+1), 2.], [complex(2**53+1), 2+0j], [2**53+1, 2]),
):
    parameters["value"] = value
    assert native.vector(parameters) == real
    assert native.complex_vector(parameters) == complex_values
    assert native.integer_vector(parameters) == integers

for value, expected in (
    (["Energy", "Stiffness"], "Energy,Stiffness"),
    ([1, "two", 3.0], "1,two,3.0"),
    ((True, False), "True,False"),
    ([None, {"key": 1}], "None,{'key': 1}"),
    (["", "middle", ""], ",middle,"),
    ([], ""),
):
    parameters["value"] = value
    assert native.text(parameters) == expected
assert native.native_text() == ",middle,"

with tempfile.TemporaryDirectory() as directory:
    with hdf5.archive(directory + "/parameters.h5", "w") as archive:
        archive["value"] = np.array([2.0, 4.0])
        archive["metadata"] = {"label": "test", "matrix": np.ones((2, 3))}
        native.load(clone, archive)
        assert native.threaded_vector(clone) == [2.0, 4.0]
        assert clone["metadata"]["matrix"].shape == (2, 3)
        clone["value"] *= 2
        assert native.vector(clone) == [4.0, 8.0]

    # A C++-created params object has no binding-owned checkpoint decoder.
    # Boolean/mixed lists use numbered HDF5 groups, which that native reader
    # must understand as well. More than ten elements catches lexical order.
    for i, values in enumerate((
        [True, False], [True, 2, 3.5], [1, "2", 3.5],
        [True] + [j + 0.5 for j in range(1, 12)],
        [True, np.int8(-7), np.int64(2**53 + 1)],
    )):
        with hdf5.archive(directory + f"/native-list-{i}.h5", "w") as archive:
            archive["parameters"] = ngs.params({"value": values})
            archive.set_context("/parameters")
            native_parameters = native.empty_vectors()
            native.load(native_parameters, archive)
            assert native.vector(native_parameters) == [float(value) for value in values]
            assert native.integer_vector(native_parameters) == [int(value) for value in values]
            assert archive.context == "/parameters"
            archive.set_context("/resaved")
            native.save(native_parameters, archive)
            restored = native.empty_vectors()
            native.load(restored, archive)
            assert native.vector(restored) == [float(value) for value in values]
            assert native.integer_vector(restored) == [int(value) for value in values]
            assert archive.context == "/resaved"
            restored["value"][0] = 7
            assert native.vector(restored)[0] == 7

    with hdf5.archive(directory + "/native-complex-list.h5", "w") as archive:
        archive["parameters"] = ngs.params({"value": [True, 2+3j, 4.]})
        archive.set_context("/parameters")
        native_parameters = native.empty_vectors()
        native.load(native_parameters, archive)
        assert native.complex_vector(native_parameters) == [1+0j, 2+3j, 4+0j]
        native.save(native_parameters, archive)
        native.load(native_parameters, archive)
        assert native.complex_vector(native_parameters) == [1+0j, 2+3j, 4+0j]

# Native destruction on a thread that started without the GIL must release
# the Python value safely, including when it owns the final reference.
value = np.ones(3)
reference = weakref.ref(value)
parameters = ngs.params({"value": value})
del value
native.destroy_on_worker(parameters)
gc.collect()
assert reference() is None
print("native parameter contracts: ok")
