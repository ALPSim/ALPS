"""Check dtype-dependent user operations, not just equal numeric values."""

import numpy as np
import pytest

from pyalps import hdf5, ngs


@pytest.mark.parametrize("dtype", [np.bool_, np.int8])
@pytest.mark.parametrize("shape", [(), (3,), (2, 3), (2, 1, 3), (0,), (2, 0)])
@pytest.mark.parametrize("attribute", [False, True])
def test_boolean_and_signed_byte_round_trip(tmp_path, dtype, shape, attribute):
    size = int(np.prod(shape))
    value = np.arange(size, dtype=np.int8).reshape(shape).astype(dtype)
    if value.ndim > 1:
        value = np.asfortranarray(value)
    value.flags.writeable = False
    path = "/group/@value" if attribute else "/value"
    filename = str(tmp_path / "dtype.h5")
    with hdf5.archive(filename, "w") as archive:
        archive.create_group("/group")
        archive[path] = value
    with hdf5.archive(filename, "r") as archive:
        actual = archive[path]
        if shape:
            assert actual.dtype == dtype
            assert actual.shape == shape
        else:
            assert type(actual) is (bool if dtype == np.bool_ else int)
        np.testing.assert_array_equal(actual, value)


def test_boolean_mask_remains_a_mask_after_reload(tmp_path):
    with hdf5.archive(str(tmp_path / "mask.h5"), "w") as archive:
        archive["mask"] = np.array([True, False, True])
        np.testing.assert_array_equal(np.array([10, 20, 30])[archive["mask"]], [10, 30])
        # Cover the native vector<bool> writer as well as ndarray dispatch.
        archive["parameters"] = ngs.params({"mask": [True, False, True]})
        np.testing.assert_array_equal(
            np.array([10, 20, 30])[archive["parameters/mask"]], [10, 30]
        )


@pytest.mark.parametrize("attribute", [False, True])
def test_overwriting_same_storage_type_updates_dtype_marker(tmp_path, attribute):
    path = "/group/@value" if attribute else "/value"
    with hdf5.archive(str(tmp_path / "overwrite.h5"), "w") as archive:
        archive.create_group("/group")
        for dtype in (np.bool_, np.int8, np.bool_, np.int8):
            archive[path] = np.array([0, 1], dtype=dtype)
            assert archive[path].dtype == dtype


def test_unmarked_legacy_boolean_dataset(tmp_path):
    with hdf5.archive(str(tmp_path / "legacy.h5"), "w") as archive:
        archive["mask"] = np.array([True, False, True])
        archive.delete_attribute("mask/@__alps_type__")
        mask = archive["mask"]
        assert mask.dtype == np.bool_
        np.testing.assert_array_equal(np.arange(3)[mask], [0, 2])
