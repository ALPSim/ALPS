 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 #                                                                                 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   #
 #               2016 - 2016 by Michele Dolfi <dolfim@phys.ethz.ch>                #
 #                                                                                 #
 # ALPS Project: https://alps.comp-phys.org/                                       #
 # SPDX-License-Identifier: MIT                                                    #
 #                                                                                 #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Assertion-based HDF5 round-trip test. The expectations encode the
# legacy Boost.Python on-disk behaviour recorded in pyhdf5io.output:
# exact-type homogeneous lists keep their element type on disk
# ([1, 2, 3] stays int32), bool/mixed/ragged lists become groups that
# read back as lists, and equal-shape numpy-array lists stack into one
# dataset with numpy's dtype.

import os
import tempfile

import numpy as np
import pyalps.hdf5 as hdf5


def _write_all(ar):
    a = np.array([1, 2, 3])
    b = np.array([1.1, 2.0, 3.5])
    c = np.array([1.1 + 1j, 2.0j, 3.5])
    d = {"a": a, 2 + 3j: "foo"}

    ar["/list"] = [1, 2, 3]
    ar["/list2"] = [[[1, 2], [3, 4]], [[1, 2], [3, 4]], [[1, 2], [3, 4]], [[1, 2], [3, 4]]]
    ar["/tuple"] = (1, 2, 3)
    ar["/dict"] = {"scalar": 1, "numpy": a, "numpycpx": c, "list": [1, 2, 3], "string": "str", 1: 1, 4: d}
    ar["/numpy"] = a
    ar["/numpy2"] = b
    ar["/numpy3"] = c
    ar["/numpyel"] = a[0]
    ar["/numpyel2"] = b[0]
    ar["/numpyel3"] = c[0]
    ar["/int"] = int(1)
    ar["/long"] = 1
    ar["/double"] = float(1)
    ar["/complex"] = complex(1, 1)
    ar["/string"] = "str"
    ar["/stringlist"] = ['a', 'list', 'of', 'strings']
    ar["/inhomogenious"] = [[1, 2, 3], a, "gurke", [[a, 2, 3], ["x", complex(1, 1)]]]
    ar["/inhomogenious2"] = [[[1, 2], [3, 4]], [[1, 2], [3, 4]], [[1, 2], [3, 4]], [[1, 2], [3]]]
    ar["/inhomogenious3"] = [np.arange(3), np.arange(5)]
    ar["/inhomogenious4"] = [np.arange(3), 10 * np.arange(3)]
    ar["/inhomogenious5"] = [list(range(3)), list(range(5)), list(range(3))]
    ar["/numpylist1"] = [np.arange(5), np.arange(5, 10)]
    ar["/numpylist2"] = [np.arange(5), np.arange(10)]
    # regression cases for the nanobind save path
    ar["/floatlist"] = [1.5, 2.5]
    ar["/cplxlist"] = [1 + 1j, 2j]
    ar["/boollist"] = [True, False]
    ar["/mixedlist"] = [1, 2.5]
    ar["/biglist"] = [2 ** 40, 2 ** 41]
    ar["/npscalars"] = list(np.arange(3))          # numpy.int64 scalars
    ar["/npboollist"] = list(np.array([True, False]))
    ar["/numpylist3"] = [np.arange(3), [3, 4, 5]]  # ndarray/list mix, legacy stacked
    ar["/boolmix"] = [np.arange(2), [True, False]]  # bool leaves veto stacking
    ar["/longmixed"] = [1, 2.5, "x"] + list(range(10))  # 13-child group
    ar["/emptylist"] = []
    ar["/shrink"] = [1, "a", "b"]
    ar["/shrink"] = [1, "a"]                       # group re-save must drop stale children


def _assert_int_array(value, expected, dtype=np.int32):
    assert isinstance(value, np.ndarray), repr(value)
    assert value.dtype == dtype, "expected %s, got %s" % (dtype, value.dtype)
    np.testing.assert_array_equal(value, expected)


def test_hdf5io():
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "pyngs.h5")
        ar = hdf5.archive(path, "w")
        _write_all(ar)
        del ar

        ar = hdf5.archive(path, "r")

        assert len(ar.list_children("/")) == 35

        # homogeneous int lists/tuples keep the int element type on disk
        _assert_int_array(ar["/list"], [1, 2, 3])
        _assert_int_array(ar["/tuple"], [1, 2, 3])
        list2 = ar["/list2"]
        assert list2.dtype == np.int32
        assert list2.shape == (4, 2, 2)
        np.testing.assert_array_equal(list2, [[[1, 2], [3, 4]]] * 4)

        # dict → group keyed by stringified keys
        d = ar["/dict"]
        assert sorted(d.keys()) == ["1", "4", "list", "numpy", "numpycpx", "scalar", "string"]
        assert d["scalar"] == 1 and d["1"] == 1 and d["string"] == "str"
        _assert_int_array(d["list"], [1, 2, 3])
        np.testing.assert_array_equal(d["numpy"], [1, 2, 3])
        np.testing.assert_allclose(d["numpycpx"], [1.1 + 1j, 2.0j, 3.5])
        assert d["4"]["(2+3j)"] == "foo"
        np.testing.assert_array_equal(d["4"]["a"], [1, 2, 3])

        # numpy arrays and scalars round-trip
        np.testing.assert_array_equal(ar["/numpy"], [1, 2, 3])
        assert np.issubdtype(ar["/numpy"].dtype, np.integer)
        np.testing.assert_allclose(ar["/numpy2"], [1.1, 2.0, 3.5])
        np.testing.assert_allclose(ar["/numpy3"], [1.1 + 1j, 2.0j, 3.5])
        assert ar["/numpyel"] == 1
        assert abs(ar["/numpyel2"] - 1.1) < 1e-12
        assert ar["/numpyel3"] == 1.1 + 1j

        # python scalars keep their types
        assert type(ar["/int"]) is int and ar["/int"] == 1
        assert type(ar["/long"]) is int and ar["/long"] == 1
        assert type(ar["/double"]) is float and ar["/double"] == 1.0
        assert type(ar["/complex"]) is complex and ar["/complex"] == 1 + 1j
        assert ar["/string"] == "str"
        assert ar["/stringlist"] == ['a', 'list', 'of', 'strings']

        # heterogeneous list → group, read back as a list
        i1 = ar["/inhomogenious"]
        assert isinstance(i1, list) and len(i1) == 4
        _assert_int_array(i1[0], [1, 2, 3])
        np.testing.assert_array_equal(i1[1], [1, 2, 3])
        assert i1[2] == "gurke"
        np.testing.assert_array_equal(i1[3][0][0], [1, 2, 3])
        assert i1[3][0][1] == 2 and i1[3][0][2] == 3
        assert i1[3][1] == ["x", 1 + 1j]

        # rectangular prefix + one ragged entry → group of matrices
        i2 = ar["/inhomogenious2"]
        assert isinstance(i2, list) and len(i2) == 4
        for entry in i2[:3]:
            assert entry.dtype == np.int32 and entry.shape == (2, 2)
            np.testing.assert_array_equal(entry, [[1, 2], [3, 4]])
        _assert_int_array(i2[3][0], [1, 2])
        _assert_int_array(i2[3][1], [3])

        # numpy-array lists: unequal shapes → group; equal shapes → stacked
        i3 = ar["/inhomogenious3"]
        assert isinstance(i3, list) and len(i3) == 2
        np.testing.assert_array_equal(i3[0], np.arange(3))
        np.testing.assert_array_equal(i3[1], np.arange(5))
        i4 = ar["/inhomogenious4"]
        assert isinstance(i4, np.ndarray) and i4.shape == (2, 3)
        assert np.issubdtype(i4.dtype, np.integer)
        np.testing.assert_array_equal(i4, [[0, 1, 2], [0, 10, 20]])
        i5 = ar["/inhomogenious5"]
        assert isinstance(i5, list) and len(i5) == 3
        for entry, size in zip(i5, (3, 5, 3)):
            _assert_int_array(entry, np.arange(size))
        nl1 = ar["/numpylist1"]
        assert isinstance(nl1, np.ndarray) and nl1.shape == (2, 5)
        assert np.issubdtype(nl1.dtype, np.integer)
        np.testing.assert_array_equal(nl1, [np.arange(5), np.arange(5, 10)])
        nl2 = ar["/numpylist2"]
        assert isinstance(nl2, list) and len(nl2) == 2
        np.testing.assert_array_equal(nl2[0], np.arange(5))
        np.testing.assert_array_equal(nl2[1], np.arange(10))

        # regression: homogeneous float / complex lists keep their type
        fl = ar["/floatlist"]
        assert fl.dtype == np.float64
        np.testing.assert_allclose(fl, [1.5, 2.5])
        cl = ar["/cplxlist"]
        assert cl.dtype == np.complex128
        np.testing.assert_allclose(cl, [1 + 1j, 2j])

        # regression: bool and mixed-type lists follow the legacy
        # per-element group behaviour instead of silently widening
        assert ar["/boollist"] == [True, False]
        ml = ar["/mixedlist"]
        assert ml == [1, 2.5]
        assert type(ml[0]) is int and type(ml[1]) is float

        # regression: out-of-int32-range values widen to int64, not float
        bl = ar["/biglist"]
        assert np.issubdtype(bl.dtype, np.integer)
        np.testing.assert_array_equal(bl, [2 ** 40, 2 ** 41])

        # regression: numpy-scalar lists keep their dtype in one dataset
        # (numpy.int64 etc. were vectorizable in the legacy build)
        nps = ar["/npscalars"]
        assert isinstance(nps, np.ndarray) and np.issubdtype(nps.dtype, np.integer)
        np.testing.assert_array_equal(nps, [0, 1, 2])
        # HDF5 has no native bool: bool arrays are stored (and read
        # back) as their int8 storage type; only the values survive
        npb = ar["/npboollist"]
        assert isinstance(npb, np.ndarray)
        assert npb.dtype == np.bool_ or npb.dtype == np.int8
        np.testing.assert_array_equal(npb, [1, 0])

        # regression: rectangular ndarray/list mixes stack, like legacy
        nl3 = ar["/numpylist3"]
        assert isinstance(nl3, np.ndarray) and nl3.shape == (2, 3)
        np.testing.assert_array_equal(nl3, [[0, 1, 2], [3, 4, 5]])

        # regression: plain bools among the leaves veto stacking — numpy
        # would silently promote True to 1
        bm = ar["/boolmix"]
        assert isinstance(bm, list) and len(bm) == 2
        np.testing.assert_array_equal(bm[0], [0, 1])
        assert bm[1] == [True, False]

        # regression: a group-saved list with more than ten elements
        # still loads as a list, in order (children come back from HDF5
        # lexicographically)
        lm = ar["/longmixed"]
        assert lm == [1, 2.5, "x"] + list(range(10)), lm

        # regression: empty lists stay integer-typed datasets
        el = ar["/emptylist"]
        assert len(el) == 0

        # regression: re-saving a group-shaped list drops stale children
        assert ar["/shrink"] == [1, "a"]

        del ar


def test_hdf5_empty_dict_roundtrip():
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "empty-dict.h5")
        ar = hdf5.archive(path, "w")
        ar["/value"] = {}
        del ar

        ar = hdf5.archive(path, "r")
        value = ar["/value"]
        assert type(value) is dict
        assert value == {}
        del ar


def test_hdf5_dict_key_roundtrip():
    expected = {
        "a/b": 1,
        "a": {"b": 2},
        "amp&key": 3,
        "entity&#47;": 4,
    }
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "dict-keys.h5")
        ar = hdf5.archive(path, "w")
        ar["/value"] = expected
        ar.create_group("/rawamp")
        ar["/rawamp/literal&child"] = 5
        del ar

        ar = hdf5.archive(path, "r")
        assert ar["/value"] == expected
        # A raw ampersand from a non-pyalps HDF5 producer is not an
        # encoded path entity and must remain literal.
        assert ar["/rawamp"] == {"literal&child": 5}
        del ar


def test_hdf5_nested_numpy_scalar_vectorization():
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "nested-numpy-scalars.h5")
        ar = hdf5.archive(path, "w")
        ar["/rectangular"] = [
            [np.int16(1), np.int16(2)],
            [np.int16(3), np.int16(4)],
        ]
        ar["/ragged"] = [
            [np.int16(1)],
            [np.int16(2), np.int16(3)],
        ]
        ar["/mixed"] = [
            [np.int16(1), np.int16(2)],
            [np.int32(3), np.int32(4)],
        ]
        ar["/arrayandscalar"] = [np.arange(2), np.int64(3)]
        del ar

        ar = hdf5.archive(path, "r")
        rectangular = ar["/rectangular"]
        assert isinstance(rectangular, np.ndarray)
        assert rectangular.dtype == np.int16
        assert rectangular.shape == (2, 2)
        np.testing.assert_array_equal(rectangular, [[1, 2], [3, 4]])

        # The new recursive case must not widen its acceptance: ragged
        # trees, mixed NumPy scalar dtypes, and sequence/scalar mixtures
        # retain the existing group representation.
        ragged = ar["/ragged"]
        assert isinstance(ragged, list) and len(ragged) == 2
        assert all(row.dtype == np.int16 for row in ragged)
        np.testing.assert_array_equal(ragged[0], [1])
        np.testing.assert_array_equal(ragged[1], [2, 3])
        mixed = ar["/mixed"]
        assert isinstance(mixed, list) and len(mixed) == 2
        assert mixed[0].dtype == np.int16
        assert mixed[1].dtype == np.int32
        array_and_scalar = ar["/arrayandscalar"]
        assert isinstance(array_and_scalar, list)
        np.testing.assert_array_equal(array_and_scalar[0], [0, 1])
        assert array_and_scalar[1] == 3
        del ar


def test_hdf5_zero_dimensional_and_zero_extent_arrays():
    scalar_cases = [
        np.array(True),
        np.array(-3, dtype=np.int32),
        np.array(2**40, dtype=np.int64),
        np.array(1.25, dtype=np.float32),
        np.array(1.25, dtype=np.float64),
        np.array(1 + 2j, dtype=np.complex64),
        np.array(1 + 2j, dtype=np.complex128),
    ]
    empty_cases = [
        np.empty((0,), dtype=np.float64),
        np.empty((0, 2), dtype=np.int32),
        np.empty((2, 0), dtype=np.int32),
        np.empty((2, 0, 3), dtype=np.complex128),
    ]

    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "zero-shapes.h5")
        with hdf5.archive(path, "w") as ar:
            for index, value in enumerate(scalar_cases):
                ar[f"/scalar/{index}"] = value
            for index, value in enumerate(empty_cases):
                ar[f"/empty/{index}"] = value

        with hdf5.archive(path, "r") as ar:
            for index, expected in enumerate(scalar_cases):
                actual = ar[f"/scalar/{index}"]
                assert np.asarray(actual).shape == ()
                assert actual == expected.item()
                if expected.dtype == np.complex64:
                    assert np.asarray(actual).dtype == np.complex64
            for index, expected in enumerate(empty_cases):
                actual = ar[f"/empty/{index}"]
                assert isinstance(actual, np.ndarray)
                assert actual.shape == expected.shape
                assert actual.dtype == expected.dtype


def test_hdf5_complex_array_precision_roundtrip():
    values = [
        np.array([1 + 2j, 3 + 4j], dtype=np.complex64),
        np.array([[1 + 2j], [3 + 4j]], dtype=np.complex64),
        np.array([1 + 2j, 3 + 4j], dtype=np.complex128),
    ]
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "complex-precision.h5")
        with hdf5.archive(path, "w") as ar:
            for index, value in enumerate(values):
                ar[f"/{index}"] = value
        with hdf5.archive(path, "r") as ar:
            for index, expected in enumerate(values):
                actual = ar[f"/{index}"]
                assert actual.dtype == expected.dtype
                np.testing.assert_array_equal(actual, expected)


def test_hdf5_strided_array_roundtrip():
    base = np.arange(24, dtype=np.float64).reshape(4, 6)
    values = [
        base[:, ::2],
        base.T,
        base[::-1, ::-2],
        np.ma.array(base[:, ::2], mask=False),
        (base.astype(np.complex64) * (1 + 2j))[::2, 1::2],
    ]
    assert all(not value.flags.c_contiguous for value in values)
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "strided.h5")
        with hdf5.archive(path, "w") as ar:
            for index, value in enumerate(values):
                ar[f"/{index}"] = value
        with hdf5.archive(path, "r") as ar:
            for index, expected in enumerate(values):
                actual = ar[f"/{index}"]
                assert actual.dtype == expected.dtype
                np.testing.assert_array_equal(actual, expected)


def test_hdf5_non_native_array_error_is_actionable():
    value = np.arange(4, dtype=np.int32).byteswap().view(np.dtype(">i4"))
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "non-native.h5")
        with hdf5.archive(path, "w") as ar:
            try:
                ar["/value"] = value
                raise AssertionError("non-native arrays must be rejected")
            except RuntimeError as error:
                assert "not native" in str(error)


if __name__ == "__main__":
    test_hdf5io()
    test_hdf5_empty_dict_roundtrip()
    test_hdf5_dict_key_roundtrip()
    test_hdf5_nested_numpy_scalar_vectorization()
    test_hdf5_zero_dimensional_and_zero_extent_arrays()
    test_hdf5_complex_array_precision_roundtrip()
    test_hdf5_strided_array_roundtrip()
    test_hdf5_non_native_array_error_is_actionable()
    print("SUCCESS")
