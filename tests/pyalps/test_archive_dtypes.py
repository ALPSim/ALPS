"""Check dtype-dependent user operations, not just equal numeric values."""

import numpy as np
import pytest

from pyalps import hdf5, ngs


@pytest.mark.parametrize("depth", [2, 3])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize(
    "second_row,dtype",
    [
        ([np.int32(3), np.int32(4)], np.int32),
        ([np.float64(3), np.float64(4)], np.float64),
        ([np.complex128(3), np.complex128(4)], np.complex128),
        (np.array([3, 4], dtype=np.int32), np.int32),
    ],
)
def test_rectangular_python_numpy_rows_keep_array_contract(
    tmp_path, second_row, dtype, reverse, depth
):
    # Both rows have the same ALPS storage dtype. Plain Python integers use
    # int32, even when NumPy's platform-default integer is int64.
    scalar = {np.int32: int, np.float64: float, np.complex128: complex}[dtype]
    value = [[scalar(1), scalar(2)], second_row]
    if reverse:
        value.reverse()
    if depth == 3:
        value = [value, value]
    expected = np.asarray(value, dtype=dtype)
    filename = str(tmp_path / "rectangular.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["table"] = value
    with hdf5.archive(filename, "r") as archive:
        restored = archive["table"]
        assert archive.is_data("table")
    # Value-only NumPy comparisons coerce lists and miss changed arithmetic.
    assert isinstance(restored, np.ndarray)
    assert restored.dtype == expected.dtype and restored.shape == expected.shape
    np.testing.assert_array_equal(restored[..., 0], expected[..., 0])
    np.testing.assert_array_equal(restored * 2, expected * 2)


@pytest.mark.parametrize(
    "value",
    [
        [[1, np.int32(2)], [np.int32(3), np.int32(4)]],
        [[1, 2], [np.int32(3)]],
        [[True, False], [np.int32(3), np.int32(4)]],
        [[2**53 + 1, 2**53 + 3], [np.int32(3), np.int32(4)]],
    ],
)
def test_incompatible_python_numpy_rows_keep_groups(tmp_path, value):
    with hdf5.archive(str(tmp_path / "groups.h5"), "w") as archive:
        archive["table"] = value
        assert archive.is_group("table")
        restored = archive["table"]
    assert isinstance(restored, list)
    for actual, expected in zip(restored, value):
        assert list(actual) == list(expected)


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("nested", [False, True])
@pytest.mark.parametrize(
    "other_dtype,integer_dtype,integers",
    [
        (np.int64, np.uint64, [2**63 + 1, 2**63 + 3]),
        (np.uint64, np.int64, [2**53 + 1, 2**53 + 3]),
        (np.float32, np.int64, [2**53 + 1, 2**53 + 3]),
        (np.float64, np.int64, [-(2**53) - 1, -(2**53) - 3]),
        (np.complex64, np.uint64, [2**64 - 1, 2**64 - 3]),
        (np.complex128, np.int64, [2**63 - 1, -(2**63) + 1]),
        (None, np.int64, [2**53 + 1, 2**53 + 3]),
    ],
)
def test_lossy_mixed_rows_preserve_exact_integers(
    tmp_path, other_dtype, integer_dtype, integers, reverse, nested
):
    other = np.array([3, 4], dtype=other_dtype) if other_dtype else [3.0, 4.0]
    # Exercise a noncontiguous, read-only integer row too.
    row = np.repeat(np.array(integers, dtype=integer_dtype), 2)[::2]
    row.flags.writeable = False
    rows = [other, row]
    if reverse:
        rows.reverse()
    value = [rows, rows] if nested else rows
    filename = str(tmp_path / "exact-integers.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["table"] = np.zeros((2, 2))  # Replace an existing dataset.
        archive["table"] = value
    with hdf5.archive(filename, "r") as archive:
        assert archive.is_group("table")
        restored = archive["table"]
    for actual in restored if nested else [restored]:
        integer_row = actual[0 if reverse else 1]
        assert integer_row.dtype == integer_dtype
        # NumPy mixed int/float equality would silently round the expected
        # integers too. Compare independently converted Python ints.
        assert [int(x) for x in integer_row] == integers
        np.testing.assert_array_equal(actual[1 if reverse else 0], other)


@pytest.mark.parametrize("scalar", [int, np.int64])
@pytest.mark.parametrize("nested", [False, True])
def test_integer_scalar_rows_are_checked_before_stacking(tmp_path, scalar, nested):
    integers = [2**53 + 1, 2**53 + 3]
    rows = [np.array([1.5, 2.5]), tuple(scalar(x) for x in integers)]
    value = [rows] if nested else rows
    filename = str(tmp_path / "scalar-rows.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["table"] = value
    with hdf5.archive(filename, "r") as archive:
        restored = archive["table"]
    if nested:
        restored = restored[0]
    assert [int(x) for x in restored[1]] == integers


@pytest.mark.parametrize(
    "other_dtype", [np.uint64, np.float32, np.float64, np.complex64]
)
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("integers", [[0, 2**40], [2**53, 2**53 + 2], [-(2**63), 0]])
def test_lossless_mixed_rows_still_stack(tmp_path, other_dtype, reverse, integers):
    other = np.array([3, 4], dtype=other_dtype)
    rows = [other, np.array(integers, dtype=np.int64)]
    if reverse:
        rows.reverse()
    filename = str(tmp_path / "lossless.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["table"] = {"old": 1}  # Replace an existing group.
        archive["table"] = rows
    with hdf5.archive(filename, "r") as archive:
        assert archive.is_data("table")
        restored = archive["table"]
    assert isinstance(restored, np.ndarray)
    assert restored.shape == (2, 2)
    assert restored.dtype == np.asarray(rows).dtype
    assert [int(x.real) for x in restored[0 if reverse else 1]] == integers


def test_lossless_mixed_rows_allow_nan_and_infinity(tmp_path):
    with hdf5.archive(str(tmp_path / "nonfinite.h5"), "w") as archive:
        archive["table"] = [np.array([np.nan, np.inf]), np.array([2**53, 0])]
        assert archive.is_data("table")
        restored = archive["table"]
    assert np.isnan(restored[0, 0]) and np.isposinf(restored[0, 1])
    assert [int(x) for x in restored[1]] == [2**53, 0]


@pytest.mark.parametrize("shape", [(), (0,), (2, 0)])
def test_mixed_integer_array_shapes(tmp_path, shape):
    rows = [np.zeros(shape, dtype=np.int64), np.full(shape, 2**64 - 1, dtype=np.uint64)]
    filename = str(tmp_path / "shapes.h5")
    with hdf5.archive(filename, "w") as archive:
        archive["table"] = rows
    with hdf5.archive(filename, "r") as archive:
        restored = archive["table"]
    if shape == ():
        assert restored == [0, 2**64 - 1]
    else:
        assert isinstance(restored, np.ndarray)
        assert restored.shape == (2, *shape)


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
