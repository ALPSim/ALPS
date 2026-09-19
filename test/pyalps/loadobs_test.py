# ****************************************************************************
#
# ALPS Project: Algorithms and Libraries for Physics Simulations
#
# ALPS Libraries
#
# Copyright (C) 2010 by Lukas Gamper <gamperl@gmail.com>
#                      Matthias Troyer <troyer@itp.phys.ethz.ch>
#               2026 by the ALPS collaboration
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
#
# ****************************************************************************

# Formerly loadobs.py (plus loadobs.cpp, its C++ twin). Both *read* a
# loadobs.h5 that neither of them created and that no longer ships anywhere in
# the tree, so both had been dead for years -- the Python half was additionally
# uncollectable by pytest, and the C++ half was orphaned when
# test/pyalps/CMakeLists.txt was removed with the in-tree bindings build.
#
# The surface they meant to cover is real: walk /simulation/results in an
# archive, decide scalar vs vector from the shape of mean/value, and load each
# entry into the matching mcdata type by path. This version builds its own
# fixture from pyalps.alea, so it needs no external file and no C++ helper.
#
# It also pins down the encoding rule both originals got wrong. list_children
# returns each name exactly as stored -- already HDF5-encoded -- so it belongs
# in a path verbatim, and hdf5_name_decode is what turns it back into the
# observable's real name. Both old versions instead re-encoded what
# list_children handed them (Python via hdf5_name_encode, C++ via
# encode_segment), which double-encodes any name containing '/' or '&':
# "Energy Density / Site" became "Energy Density &#38;#47; Site" and the
# lookup missed. Their fixture evidently held only simple names.

import numpy as np
import pytest

import pyalps
import pyalps.alea as alea
import pyalps.hdf5 as hdf5


RESULTS = "/simulation/results"

SCALARS = {
    "Energy": (-1.2345, 0.0067),
    # '^' and '/' are exactly why hdf5_name_encode exists: they would
    # otherwise be read as HDF5 path syntax.
    "Magnetization^2": (0.3721, 0.0014),
    "Energy Density / Site": (-0.6172, 0.0033),
}

VECTORS = {
    "Correlations": ([1.0, 0.5, 0.25], [0.01, 0.02, 0.03]),
    "Green's Function": ([0.9, -0.4], [0.05, 0.06]),
}


@pytest.fixture
def results_archive(tmp_path):
    """An archive holding scalar and vector observables under /simulation/results."""
    path = str(tmp_path / "loadobs.h5")

    for name, (mean, error) in SCALARS.items():
        observable = alea.MCScalarData(mean, error)
        observable.save(path, RESULTS + "/" + pyalps.hdf5_name_encode(name))

    for name, (mean, error) in VECTORS.items():
        observable = alea.MCVectorData(np.array(mean), np.array(error))
        observable.save(path, RESULTS + "/" + pyalps.hdf5_name_encode(name))

    return path


def test_list_children_returns_stored_encoded_names(results_archive):
    """Names come back encoded; decoding is what recovers the real name."""
    with hdf5.archive(results_archive, "r") as archive:
        children = archive.list_children(RESULTS)

    expected = list(SCALARS) + list(VECTORS)
    assert sorted(children) == sorted(pyalps.hdf5_name_encode(n) for n in expected)
    assert sorted(pyalps.hdf5_name_decode(c) for c in children) == sorted(expected)

    # Re-encoding a name from list_children is the bug the originals had.
    encoded_twice = pyalps.hdf5_name_encode("Energy Density &#47; Site")
    assert encoded_twice != "Energy Density &#47; Site"
    with hdf5.archive(results_archive, "r") as archive:
        assert not archive.is_group(RESULTS + "/" + encoded_twice)


def test_name_encoding_roundtrips_through_the_archive(results_archive):
    for name in list(SCALARS) + list(VECTORS):
        encoded = pyalps.hdf5_name_encode(name)
        assert pyalps.hdf5_name_decode(encoded) == name
        with hdf5.archive(results_archive, "r") as archive:
            assert archive.is_group(RESULTS + "/" + encoded)


def test_scalar_and_vector_observables_are_distinguishable(results_archive):
    """The scalar/vector decision the original test drove its dispatch from."""
    with hdf5.archive(results_archive, "r") as archive:
        for name in SCALARS:
            path = RESULTS + "/" + pyalps.hdf5_name_encode(name) + "/mean/value"
            assert archive.is_scalar(path), name
        for name in VECTORS:
            path = RESULTS + "/" + pyalps.hdf5_name_encode(name) + "/mean/value"
            assert not archive.is_scalar(path), name


def test_load_dispatches_on_shape_and_restores_the_values(results_archive):
    """The whole original loop, now asserting instead of printing."""
    with hdf5.archive(results_archive, "r") as archive:
        loaded = {}
        for stored in archive.list_children(RESULTS):
            # `stored` is already encoded: use it verbatim in the path, and
            # decode only to recover the observable's name.
            path = RESULTS + "/" + stored
            if archive.is_scalar(path + "/mean/value"):
                observable = alea.MCScalarData()
            else:
                observable = alea.MCVectorData()
            observable.load(results_archive, path)
            loaded[pyalps.hdf5_name_decode(stored)] = observable

    assert sorted(loaded) == sorted(list(SCALARS) + list(VECTORS))

    for name, (mean, error) in SCALARS.items():
        assert isinstance(loaded[name], alea.MCScalarData), name
        assert loaded[name].mean == pytest.approx(mean), name
        assert loaded[name].error == pytest.approx(error), name

    for name, (mean, error) in VECTORS.items():
        assert isinstance(loaded[name], alea.MCVectorData), name
        np.testing.assert_allclose(loaded[name].mean, mean)
        np.testing.assert_allclose(loaded[name].error, error)


def test_loaded_observables_are_printable(results_archive):
    observable = alea.MCScalarData()
    observable.load(results_archive, RESULTS + "/" + pyalps.hdf5_name_encode("Energy"))
    assert "+/-" in str(observable)
