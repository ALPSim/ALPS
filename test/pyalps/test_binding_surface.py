#!/usr/bin/env python3
# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Lock the public pyalps extension surface after the nanobind migration."""

from __future__ import annotations

import copy
import importlib
import os
import tempfile
from types import SimpleNamespace

import numpy as np


def test_extension_import_surface():
    import pyalps
    import pyalps.cxx as cxx

    expected = {
        "pyalea_c",
        "pymcdata_c",
        "pytools_c",
        "pyngsparams_c",
        "pyngshdf5_c",
        "pyngsbase_c",
        "pyngsobservable_c",
        "pyngsobservables_c",
        "pyngsresult_c",
        "pyngsresults_c",
        "pyngsapi_c",
        "pyngsrandom01_c",
        "pyngsaccumulator_c",
    }
    assert pyalps is not None
    assert expected <= set(vars(cxx))


def test_cross_module_parameter_archive_and_rng_roundtrip():
    from pyalps.cxx import pyngshdf5_c, pyngsparams_c, pyngsrandom01_c

    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "surface.h5")
        params = pyngsparams_c.params()
        params["integer"] = 42
        params["real"] = 3.25
        params["flag"] = True
        params["text"] = "nanobind-lock"

        rng = pyngsrandom01_c.random01(91)
        for _ in range(7):
            rng()

        archive = pyngshdf5_c.hdf5_archive_impl(path, "w")
        archive.create_group("/parameters")
        archive.set_context("/parameters")
        params.save(archive)
        archive.set_context("/")
        rng.save(archive)
        del archive

        loaded = pyngsparams_c.params()
        restored_rng = pyngsrandom01_c.random01(0)
        archive = pyngshdf5_c.hdf5_archive_impl(path, "r")
        archive.set_context("/parameters")
        loaded.load(archive)
        archive.set_context("/")
        restored_rng.load(archive)
        del archive

        assert sorted(loaded) == sorted(params)
        assert int(loaded["integer"]) == 42
        assert float(loaded["real"]) == 3.25
        assert bool(loaded["flag"]) is True
        assert str(loaded["text"]) == "nanobind-lock"
        assert [rng() for _ in range(5)] == [restored_rng() for _ in range(5)]


def test_alea_numpy_and_mcdata_operators():
    from pyalps.cxx.pyalea_c import MCScalarTimeseries, RealObservable, mean, size
    from pyalps.cxx.pymcdata_c import MCScalarData

    observable = RealObservable("energy")
    for sample in (0.9, 1.0, 1.1, 1.0):
        observable << sample
    assert observable.count == 4
    assert abs(observable.mean - 1.0) < 1e-12
    assert observable.error >= 0

    series = MCScalarTimeseries(np.asarray([1.0, 2.0, 3.0]))
    assert size(series) == 3
    assert mean(series) == 2.0
    np.testing.assert_allclose(series.timeseries(), [1.0, 2.0, 3.0])

    first = MCScalarData(1.0, 0.1)
    second = MCScalarData(2.0, 0.2)
    total = first + second
    assert total.mean == 3.0
    assert total.error > 0
    duplicate = copy.deepcopy(total)
    assert duplicate.mean == total.mean
    assert duplicate.error == total.error


def test_alea_mcanalyze_surface():
    from pyalps import alea
    from pyalps.cxx.pyalea_c import (
        MCScalarTimeseries,
        MCScalarTimeseriesView,
        MCVectorTimeseries,
        StdPairDouble,
        integrated_autocorrelation_time,
        size,
    )
    from pyalps.cxx.pytools_c import rng

    generator = rng(42)
    samples = []
    state = 0.0
    for _ in range(512):
        state = 0.9 * state + 0.1 * (2 * generator() - 1)
        samples.append(state)
    series = MCScalarTimeseries(np.asarray(samples))

    correlation = alea.autocorrelation(series, _distance=16)
    assert isinstance(correlation, MCScalarTimeseries)
    assert size(correlation) == 16
    limited = alea.autocorrelation(series, _limit=0.2)
    assert size(limited) >= 1

    head = alea.cut_head(series, _distance=100)
    tail = alea.cut_tail(series, _distance=100)
    assert isinstance(head, MCScalarTimeseriesView)
    assert size(head) == 412
    assert size(tail) == 412
    assert size(alea.cut_head(correlation, _limit=0.5)) < 16

    fit = alea.exponential_autocorrelation_time(correlation, _from=1, _to=8)
    assert isinstance(fit, StdPairDouble)
    assert fit.second < 0  # decaying autocorrelation
    ranged = alea.exponential_autocorrelation_time(correlation, _max=0.8, _min=0.2)
    assert isinstance(ranged, StdPairDouble)

    tau_from_pair = integrated_autocorrelation_time(correlation, fit)
    tau_from_tuple = integrated_autocorrelation_time(correlation, (fit.first, fit.second))
    assert tau_from_pair == tau_from_tuple
    assert tau_from_pair > 0

    assert alea.error(series) > 0
    assert alea.error(series, "binning") > 0

    vector_series = MCVectorTimeseries(np.asarray([[float(i + j) for j in range(3)] for i in range(64)]))
    vector_error = alea.error(vector_series)
    assert vector_error.shape == (3,)
    assert np.all(vector_error > 0)
    vector_correlation = alea.autocorrelation(vector_series, _distance=4)
    assert vector_correlation.timeseries().shape == (4, 3)


def test_packaged_xml_stylesheets():
    import pyalps.tools

    xsl = pyalps.tools.xslPath()
    assert os.path.basename(xsl) == "ALPS.xsl"
    assert os.path.exists(xsl)
    xml_dir = os.path.dirname(xsl)
    for name in ("lattices.xml", "models.xml", "plot2mpl.xsl"):
        assert os.path.exists(os.path.join(xml_dir, name))


def test_ngs_observable_containers():
    from pyalps import ngs

    observables = ngs.observables()
    observables.createRealObservable("magnetization")
    observables["magnetization"] << 1.5
    assert "magnetization" in observables
    assert ngs.observable2result(observables["magnetization"]).count == 1


def test_name_encoding_roundtrip():
    from pyalps.cxx.pytools_c import hdf5_name_decode, hdf5_name_encode

    for value in ("plain", "with space", "slash/inside", "café"):
        assert hdf5_name_decode(hdf5_name_encode(value)) == value


def test_accumulator_surface():
    from pyalps.cxx.pyngsaccumulator_c import error_accumulator

    accumulator = error_accumulator()
    for sample in (1.0, 2.0, 3.0):
        accumulator(sample)
    result = accumulator.result()
    assert result.count() == 3
    assert result.mean() == 2.0
    assert result.error() >= 0


def test_optional_application_extension_surface():
    for name in ("maxent_c", "dwa_c", "cthyb", "ctint"):
        module = importlib.import_module("pyalps._ext." + name)
        assert module.__name__.endswith(name)


def test_current_python_numpy_and_scipy_compatibility(monkeypatch):
    import pyalps
    import pyalps.dwa as dwa

    assert callable(dwa.thermalized)

    parsed = pyalps.stringListToList("[1,[2,3],4]")
    assert parsed == [[1.0], [2.0, 3.0], [4.0]]

    shared = pyalps.dict_intersect([
        {"array": np.array([1, 2]), "scalar": 3},
        {"array": np.array([1, 2]), "scalar": 3},
    ])
    np.testing.assert_array_equal(shared["array"], [1, 2])
    assert shared["scalar"] == 3

    monkeypatch.setattr(
        pyalps,
        "loadTimeSeries",
        lambda *_args: np.array([1.0, 1.1, 0.9, 1.0]),
    )
    steady = pyalps.checkSteadyState(outfile="unused.h5", observable="energy")
    assert isinstance(steady["value"], (bool, np.bool_))


def test_python3_property_comparison(monkeypatch):
    import pyalps
    import pyalps.apptest as apptest

    properties = {
        "test.h5": {"vector": np.array([1, 2]), "value": 3},
        "reference.h5": {"vector": np.array([1, 2]), "value": 3},
    }

    class Loader:
        def GetProperties(self, filenames):
            return [SimpleNamespace(props=properties[filenames[0]].copy())]

    monkeypatch.setattr(pyalps.load, "Hdf5Loader", Loader)
    assert apptest.checkProperties("test.h5", "reference.h5")


if __name__ == "__main__":
    for test in (
        test_extension_import_surface,
        test_cross_module_parameter_archive_and_rng_roundtrip,
        test_alea_numpy_and_mcdata_operators,
        test_alea_mcanalyze_surface,
        test_packaged_xml_stylesheets,
        test_ngs_observable_containers,
        test_name_encoding_roundtrip,
        test_accumulator_surface,
        test_optional_application_extension_surface,
    ):
        test()
    print("pyalps binding surface: green")
