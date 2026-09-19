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

with tempfile.TemporaryDirectory() as directory:
    with hdf5.archive(directory + "/parameters.h5", "w") as archive:
        archive["value"] = np.array([2.0, 4.0])
        archive["metadata"] = {"label": "test", "matrix": np.ones((2, 3))}
        native.load(clone, archive)
        assert native.threaded_vector(clone) == [2.0, 4.0]
        assert clone["metadata"]["matrix"].shape == (2, 3)
        clone["value"] *= 2
        assert native.vector(clone) == [4.0, 8.0]

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
