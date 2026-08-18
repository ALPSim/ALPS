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
    c = np.array([1.1 + 1j, 2.0j, 3.5])
    d = {"a": a, 2 + 3j: "foo"}

    ar["/list"] = [1, 2, 3]
    ar["/list2"] = [[[1, 2], [3, 4]], [[1, 2], [3, 4]], [[1, 2], [3, 4]], [[1, 2], [3, 4]]]
    ar["/tuple"] = (1, 2, 3)
    ar["/dict"] = {"scalar": 1, "numpy": a, "numpycpx": c, "list": [1, 2, 3], "string": "str", 1: 1, 4: d}
    ar["/numpy"] = a
    ar["/numpy2"] = np.array([1.1, 2.0, 3.5])
    ar["/numpy3"] = c
    ar["/numpyel"] = a[0]
    ar["/numpyel2"] = np.array([1.1, 2.0, 3.5])[0]
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

        assert len(ar.list_children("/")) == 28

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

        del ar


if __name__ == "__main__":
    test_hdf5io()
    print("SUCCESS")
