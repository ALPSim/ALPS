#!/usr/bin/env python3
# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Lock the public pyalps extension surface after the nanobind migration."""

from __future__ import annotations

import copy
import importlib
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace

import numpy as np

from pyalps.hdf5 import archive as hdf5_archive
import pytest


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
    for name in expected:
        direct = importlib.import_module("pyalps." + name)
        assert direct is getattr(cxx, name)
        assert getattr(pyalps, name) is direct


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


def test_params_from_parameter_file():
    from pyalps.cxx.pyngsparams_c import params

    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "input.parm")
        with open(path, "w") as parameter_file:
            parameter_file.write('LATTICE="chain lattice";\nL=10;\nT=2.25;\n')
        loaded = params(path)
        assert str(loaded["LATTICE"]) == "chain lattice"
        assert int(loaded["L"]) == 10
        assert float(loaded["T"]) == 2.25


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
    for name in ("maxent_c", "cthyb", "ctint"):
        module = importlib.import_module("pyalps._ext." + name)
        assert module.__name__.endswith(name)
        assert importlib.import_module("pyalps." + name) is module

    from pyalps import cthyb, ctint, maxent_c
    assert callable(maxent_c.AnalyticContinuation)
    assert callable(cthyb.solve)
    assert callable(ctint.solve)

def test_ctqmc_solvers_restore_python_signal_handlers(tmp_path, monkeypatch):
    from pyalps import cthyb, ctint
    import pyalps.hdf5 as hdf5

    monkeypatch.chdir(tmp_path)

    delta_path = tmp_path / "delta.dat"
    delta_path.write_text("".join(f"{i} -0.5 -0.5\n" for i in range(11)))
    cthyb_params = {
        "SWEEPS": 1,
        "MAX_TIME": 1,
        "THERMALIZATION": 0,
        "SEED": 0,
        "N_MEAS": 1,
        "N_HISTOGRAM_ORDERS": 4,
        "N_ORBITALS": 2,
        "U": 1.0,
        "MU": 0.5,
        "DELTA": str(delta_path),
        "N_TAU": 10,
        "BETA": 1.0,
        "TEXT_OUTPUT": 0,
        "BASENAME": str(tmp_path / "cthyb-signal"),
    }

    ctint_input = tmp_path / "ctint-input.h5"
    archive = hdf5.archive(str(ctint_input), "w")
    bare_green = np.asarray([-1j, -0.3j, -0.2j, -0.1j])
    archive["/G0_0"] = bare_green
    archive["/G0_1"] = bare_green
    del archive
    ctint_params = {
        "SWEEPS": 1,
        "MAX_TIME": 1,
        "THERMALIZATION": 0,
        "BETA": 1.0,
        "U": 1.0,
        "MU": 0.5,
        "ALPHA": 0.5,
        "N_MATSUBARA": 4,
        "N_TAU": 4,
        "INFILE": str(ctint_input),
        "BASENAME": str(tmp_path / "ctint-signal"),
    }

    calls = []

    def python_sigint_handler(signum, frame):
        calls.append((signum, frame))

    previous_handler = signal.signal(signal.SIGINT, python_sigint_handler)
    try:
        # Run each solver twice: restoration alone is not enough if ALPS' own
        # handlers are not reinstalled for the next embedded call.
        for solver, params in ((cthyb, cthyb_params), (ctint, ctint_params)):
            for _ in range(2):
                solver.solve(params)
                assert signal.getsignal(signal.SIGINT) is python_sigint_handler
                signal.raise_signal(signal.SIGINT)
                assert calls[-1][0] == signal.SIGINT
    finally:
        signal.signal(signal.SIGINT, previous_handler)

    assert len(calls) == 4


def test_maxent_restores_python_signal_handlers(tmp_path, monkeypatch):
    """MaxEnt must hand SIGINT back to Python, like cthyb and ctint do.

    Note the assertion style: signal.getsignal() is NOT a valid check here.
    ALPS installs its handler with sigaction() behind CPython's back, so
    getsignal() keeps reporting the Python handler while the OS-level
    disposition belongs to ALPS -- an unguarded run passes a getsignal()
    check and still swallows Ctrl-C, printing "Received signal 2" instead.
    Only actually raising the signal and observing whether the Python
    handler runs detects it.
    """
    maxent = pytest.importorskip("pyalps.maxent_c")

    monkeypatch.chdir(tmp_path)

    ndat = 6
    parms = {
        "BETA": 2.0, "NDAT": ndat, "NFREQ": 20, "N_ALPHA": 2,
        "ALPHA_MIN": 0.1, "ALPHA_MAX": 1.0, "MAX_IT": 2,
        "OMEGA_MAX": 4.0, "FREQUENCY_GRID": "linear", "KERNEL": "fermionic",
        "DATASPACE": "time", "TEXT_OUTPUT": 0, "VERBOSE": 0,
        "PARTICLE_HOLE_SYMMETRY": 1, "NORM": 1.0, "MAX_TIME": 1,
        "BASENAME": str(tmp_path / "maxent-signal"),
    }
    for index in range(ndat):
        parms["X_%d" % index] = -0.5
        parms["SIGMA_%d" % index] = 0.01

    calls = []

    def python_sigint_handler(signum, frame):
        calls.append(signum)

    previous_handler = signal.signal(signal.SIGINT, python_sigint_handler)
    try:
        # Twice: restoring once is not enough if ALPS' own handlers are not
        # reinstalled for the next embedded call.
        for _ in range(2):
            maxent.AnalyticContinuation(parms)
            signal.raise_signal(signal.SIGINT)
    finally:
        signal.signal(signal.SIGINT, previous_handler)

    assert calls == [signal.SIGINT, signal.SIGINT], (
        "SIGINT was not handed back to Python after AnalyticContinuation; "
        "ALPS still owns the OS-level handler"
    )


def test_mpi4py_compatibility_surface():
    pytest.importorskip("mpi4py")
    import operator
    from pyalps import ngs
    import pyalps.mpi as mpi

    assert mpi.initialized()
    assert mpi.world.rank == mpi.rank
    assert mpi.world.size == mpi.size
    assert issubclass(mpi.Exception, Exception)
    assert mpi.Communicator().rank == mpi.rank
    assert mpi.broadcast(value={"rank": mpi.rank}, root=0) == {"rank": 0}
    assert mpi.all_gather(value=mpi.rank) == tuple(range(mpi.size))
    gathered = mpi.gather(value=mpi.rank, root=0)
    if mpi.rank == 0:
        assert gathered == tuple(range(mpi.size))
    else:
        assert gathered is None
    scattered = mpi.scatter(
        values=tuple("rank-{}".format(index) for index in range(mpi.size))
        if mpi.rank == 0 else None,
        root=0,
    )
    assert scattered == "rank-{}".format(mpi.rank)
    exchanged = mpi.all_to_all(
        values=tuple((mpi.rank, destination) for destination in range(mpi.size))
    )
    assert exchanged == tuple((source, mpi.rank) for source in range(mpi.size))
    assert mpi.reduce(value=1, op=operator.add, root=0) == (
        mpi.size if mpi.rank == 0 else None
    )
    assert mpi.all_reduce(value=1, op=operator.add) == mpi.size
    assert mpi.scan(value=mpi.rank + 1, op=operator.add) == (
        (mpi.rank + 1) * (mpi.rank + 2) // 2
    )

    subgroup = mpi.world.split(color=mpi.rank % 2, key=mpi.rank)
    assert subgroup and subgroup.rank >= 0 and subgroup.size >= 1
    mpi.world.barrier()

    # Exercise actual inter-rank transport under mpiexec, while remaining a
    # valid self-send in the ordinary one-process wheel test.
    send_to = (mpi.rank + 1) % mpi.size
    receive_from = (mpi.rank - 1) % mpi.size
    ring_request = mpi.world.isend(
        send_to, tag=172, value={"source": mpi.rank, "payload": "ring"}
    )
    ring_value, ring_status = mpi.world.recv(
        receive_from, tag=172, return_status=True
    )
    ring_request.wait()
    assert ring_value == {"source": receive_from, "payload": "ring"}
    assert ring_status.source == receive_from and ring_status.tag == 172

    # Point-to-point spelling and return_status match Boost.MPI's Python API.
    request = mpi.world.isend(mpi.rank, tag=173, value="self")
    value, status = mpi.world.recv(mpi.rank, tag=173, return_status=True)
    send_status = request.wait()
    assert value == "self"
    assert status.source == mpi.rank and status.tag == 173
    assert isinstance(send_status, mpi.Status)

    send_request = mpi.world.isend(mpi.rank, tag=174, value="async")
    receive_request = mpi.world.irecv(mpi.rank, tag=174)
    assert isinstance(send_request, mpi.Request)
    assert isinstance(receive_request, mpi.RequestWithValue)
    received, receive_status = receive_request.wait()
    send_request.wait()
    assert received == "async"
    assert receive_status.source == mpi.rank and receive_status.tag == 174

    callbacks = []
    requests = mpi.RequestList([
        mpi.world.isend(mpi.rank, tag=175, value="batch"),
        mpi.world.irecv(mpi.rank, tag=175),
    ])
    mpi.wait_all(requests, lambda result, result_status: callbacks.append(
        (result, result_status)
    ))
    assert callbacks[1][0] == "batch"
    assert callbacks[1][1].source == mpi.rank

    any_send = mpi.world.isend(mpi.rank, tag=176, value="any")
    any_requests = mpi.RequestList([mpi.world.irecv(mpi.rank, tag=176)])
    any_value, any_status, any_index = mpi.wait_any(any_requests)
    any_send.wait()
    assert (any_value, any_index) == ("any", 0)
    assert any_status.source == mpi.rank

    some_callbacks = []
    some_send = mpi.world.isend(mpi.rank, tag=177, value="some")
    some_requests = mpi.RequestList([mpi.world.irecv(mpi.rank, tag=177)])
    boundary = mpi.wait_some(
        some_requests,
        lambda result, result_status: some_callbacks.append(
            (result, result_status.source)
        ),
    )
    some_send.wait()
    assert boundary == 0
    assert some_callbacks == [("some", mpi.rank)]

    poll_send = mpi.world.isend(mpi.rank, tag=178, value="request-test")
    poll_receive = mpi.world.irecv(mpi.rank, tag=178)
    deadline = time.monotonic() + 5
    poll_result = None
    while poll_result is None and time.monotonic() < deadline:
        poll_result = poll_receive.test()
    poll_send.wait()
    assert poll_result is not None
    assert poll_result[0] == "request-test"

    any_test_send = mpi.world.isend(mpi.rank, tag=179, value="test-any")
    any_test_requests = mpi.RequestList([mpi.world.irecv(mpi.rank, tag=179)])
    deadline = time.monotonic() + 5
    any_test_result = None
    while any_test_result is None and time.monotonic() < deadline:
        any_test_result = mpi.test_any(any_test_requests)
    any_test_send.wait()
    assert any_test_result is not None
    assert (any_test_result[0], any_test_result[2]) == ("test-any", 0)

    all_test_callbacks = []
    all_test_send = mpi.world.isend(mpi.rank, tag=180, value="test-all")
    all_test_requests = mpi.RequestList([mpi.world.irecv(mpi.rank, tag=180)])
    deadline = time.monotonic() + 5
    while (not mpi.test_all(
        all_test_requests,
        lambda result, result_status: all_test_callbacks.append(
            (result, result_status.source)
        ),
    ) and time.monotonic() < deadline):
        pass
    all_test_send.wait()
    assert all_test_callbacks == [("test-all", mpi.rank)]

    some_test_callbacks = []
    some_test_send = mpi.world.isend(mpi.rank, tag=181, value="test-some")
    some_test_requests = mpi.RequestList([mpi.world.irecv(mpi.rank, tag=181)])
    deadline = time.monotonic() + 5
    some_test_boundary = len(some_test_requests)
    while some_test_boundary != 0 and time.monotonic() < deadline:
        some_test_boundary = mpi.test_some(
            some_test_requests,
            lambda result, result_status: some_test_callbacks.append(
                (result, result_status.source)
            ),
        )
    some_test_send.wait()
    assert some_test_boundary == 0
    assert some_test_callbacks == [("test-some", mpi.rank)]

    probe_send = mpi.world.isend(mpi.rank, tag=182, value="probe")
    probe_status = mpi.world.probe(mpi.rank, tag=182)
    assert probe_status.source == mpi.rank and probe_status.tag == 182
    assert mpi.world.recv(mpi.rank, tag=182) == "probe"
    probe_send.wait()

    iprobe_send = mpi.world.isend(mpi.rank, tag=183, value="iprobe")
    deadline = time.monotonic() + 5
    iprobe_status = None
    while iprobe_status is None and time.monotonic() < deadline:
        iprobe_status = mpi.world.iprobe(mpi.rank, tag=183)
    assert iprobe_status is not None
    assert iprobe_status.source == mpi.rank and iprobe_status.tag == 183
    assert mpi.world.recv(mpi.rank, tag=183) == "iprobe"
    iprobe_send.wait()

    timer = mpi.Timer()
    assert timer.elapsed >= 0
    assert 0 < timer.elapsed_min < timer.elapsed_max
    assert mpi.max_tag + 1 == mpi.collectives_tag

    # The legacy mcbase constructor accepted a communicator but did not use
    # it internally. Preserve that call shape without binding Boost.MPI.
    class Simulation(ngs.mcbase):
        def update(self):
            pass

        def measure(self):
            pass

        def fraction_completed(self):
            return 1.0

    assert isinstance(Simulation({"SEED": 1}, 42, mpi.world), ngs.mcbase)


def test_mpi_finalization_ownership():
    pytest.importorskip("mpi4py")

    # Boost.MPI finalized only an environment its Python module initialized.
    # Importing pyalps.mpi after an existing mpi4py user must therefore leave
    # that user's MPI process alive when pyalps.mpi.finalize() is called.
    externally_owned = """
from mpi4py import MPI
import pyalps.mpi as mpi
assert not mpi._initialized_here
mpi.finalize()
assert MPI.Is_initialized() and not MPI.Is_finalized()
"""
    subprocess.run([sys.executable, "-c", externally_owned], check=True)

    # Conversely, a direct pyalps.mpi import owns the initialization and its
    # explicit finalize call must release it.
    pyalps_owned = """
import pyalps.mpi as mpi
assert mpi._initialized_here
mpi.finalize()
assert mpi.finalized()
"""
    subprocess.run([sys.executable, "-c", pyalps_owned], check=True)


@pytest.mark.skipif(
    os.environ.get("PYALPS_TEST_DOWNSTREAM_EXPORT") != "1",
    reason="enabled for one wheel per platform in packaging CI",
)
def test_downstream_nanobind_simulation_export(tmp_path):
    """Build and run a consumer extension against the installed ALPS SDK."""
    repository = Path(__file__).resolve().parents[2]
    tutorial = repository / "tutorials" / "ngs" / "5_export_python"
    alps_dir = repository / "_build" / "wheel-deps" / "install" / "share" / "alps"
    build = tmp_path / "export-python-build"

    assert (alps_dir / "ALPSConfig.cmake").is_file()
    subprocess.run(
        [
            "cmake", "-S", str(tutorial), "-B", str(build),
            "-DALPS_DIR={}".format(alps_dir),
            "-DPython_EXECUTABLE={}".format(sys.executable),
        ],
        check=True,
    )
    subprocess.run(
        ["cmake", "--build", str(build), "--parallel", "2"],
        check=True,
    )

    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        filter(None, (str(build), environment.get("PYTHONPATH")))
    )
    completed = subprocess.run(
        [sys.executable, str(tutorial / "smoke_test.py")],
        capture_output=True,
        env=environment,
        text=True,
    )
    assert completed.returncode == 0, (
        "downstream exporter smoke test failed\n"
        f"stdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}"
    )
    assert "downstream nanobind export: ok" in completed.stdout


def test_current_python_numpy_and_scipy_compatibility(monkeypatch):
    import pyalps

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


def test_params_mapping_equality_and_value_ladder():
    from pyalps import ngs

    # MutableMapping equality — lost by the old hasattr-guarded shim,
    # present under the Boost.Python __bases__ inheritance
    assert ngs.params({"a": 1}) == ngs.params({"a": 1})
    assert ngs.params({"a": 1}) != ngs.params({"a": 2})
    assert ngs.params({"a": 1}) == {"a": 1}
    # and, like a Mapping with __eq__, unhashable
    try:
        hash(ngs.params({}))
        raise AssertionError("params must be unhashable")
    except TypeError:
        pass

    p = ngs.params({})
    # None is rejected with a message that says so
    try:
        p["x"] = None
        raise AssertionError("None must be rejected")
    except TypeError as error:
        assert "None" in str(error)
    # oversized integers raise instead of truncating silently —
    # inside lists too, where the double-widening fallback would
    # otherwise corrupt values beyond 2**53
    try:
        p["n"] = 2 ** 40
        raise AssertionError("2**40 must be rejected")
    except TypeError as error:
        assert "32-bit" in str(error)
    try:
        p["nl"] = [2 ** 53 + 1]
        raise AssertionError("[2**53+1] must be rejected")
    except TypeError as error:
        assert "32-bit" in str(error)
    # Homogeneous bool sequences have a native C++ representation and
    # round-trip without falling back to stored Python objects.
    p["flags"] = [True, False]
    assert p["flags"] == [True, False]
    p["npflags"] = np.array([True, False], dtype=np.bool_)
    assert p["npflags"] == [True, False]
    try:
        p["mixedflags"] = [True, 1]
        raise AssertionError("mixed bool/numeric sequences must be rejected")
    except TypeError as error:
        assert "cannot be mixed" in str(error)
    # numpy integer scalars are accepted like numpy floats are —
    # as scalars and inside lists, with the same 32-bit range policy
    p["npint"] = np.int64(8)
    assert p["npint"] == 8 and type(p["npint"]) is int
    p["npbool"] = np.bool_(True)
    assert p["npbool"] is True
    p["npfloat32"] = np.float32(1.25)
    assert p["npfloat32"] == 1.25
    p["nplongdouble"] = np.longdouble("1.125")
    assert p["nplongdouble"] == 1.125
    p["npcomplex64"] = np.complex64(1 + 2j)
    assert p["npcomplex64"] == 1 + 2j
    p["npclongdouble"] = np.clongdouble(3 + 4j)
    assert p["npclongdouble"] == 3 + 4j
    p["npbytes"] = np.bytes_(b"native")
    assert p["npbytes"] == "native"
    p["npints"] = [np.int64(1), np.int64(2)]
    assert p["npints"] == [1, 2]
    assert all(type(v) is int for v in p["npints"])
    p["nparray"] = np.array([1, 2], dtype=np.int64)
    assert p["nparray"] == [1, 2]
    p["npsubclass"] = np.ma.array([1, 2], mask=False)
    assert p["npsubclass"] == [1, 2]
    p["npfloats"] = np.array([1.5, 2.5], dtype=np.float32)
    assert p["npfloats"] == [1.5, 2.5]
    p["npcomplex"] = np.array([1 + 2j, 3 + 4j], dtype=np.complex64)
    assert p["npcomplex"] == [1 + 2j, 3 + 4j]
    p["npextendedcomplex"] = np.array(
        [1 + 2j, 3 + 4j], dtype=np.clongdouble
    )
    assert p["npextendedcomplex"] == [1 + 2j, 3 + 4j]
    p["npcomplexlist"] = [np.complex64(5 + 6j), np.clongdouble(7 + 8j)]
    assert p["npcomplexlist"] == [5 + 6j, 7 + 8j]
    p["npstrings"] = np.array(["a", "b"])
    assert p["npstrings"] == ["a", "b"]
    p["npbytestrings"] = np.array([b"a", b"b"], dtype="S1")
    assert p["npbytestrings"] == ["a", "b"]
    p["np0d"] = np.array(7, dtype=np.int64)
    assert p["np0d"] == 7
    p["emptyflags"] = np.array([], dtype=np.bool_)
    assert p["emptyflags"] == []
    try:
        p["matrix"] = np.ones((2, 2))
        raise AssertionError("multidimensional parameter arrays must be rejected")
    except TypeError as error:
        assert "multidimensional" in str(error)
    try:
        p["npbig"] = [np.int64(2 ** 40)]
        raise AssertionError("[np.int64(2**40)] must be rejected")
    except TypeError as error:
        assert "32-bit" in str(error)
    # exact-type lists round-trip with their element type
    p["ilist"] = [1, 2, 3]
    assert p["ilist"] == [1, 2, 3]
    assert all(type(v) is int for v in p["ilist"])
    p["flist"] = [1.5, 2.5]
    assert p["flist"] == [1.5, 2.5]
    p["slist"] = ["a", "b"]
    assert p["slist"] == ["a", "b"]
    # mixed numeric lists widen to double; complex scalars are stored
    p["mixed"] = [1, 2.5]
    assert p["mixed"] == [1.0, 2.5]
    p["cplx"] = 1 + 2j
    assert p["cplx"] == 1 + 2j

    # Unsupported object graphs stay unsupported: params owns only native
    # C++ values and must never keep arbitrary Python objects alive.
    for unsupported in ({"nested": 1}, object()):
        try:
            p["object"] = unsupported
            raise AssertionError("arbitrary Python objects must be rejected")
        except TypeError:
            pass


def test_params_mapping_mixins_handle_none_getitem():
    """get/pop/setdefault must honour their contracts on params.

    params.__getitem__ returns None for an undefined key rather than raising
    KeyError, so MutableMapping's mixins -- which are written against the
    KeyError contract -- silently misbehaved: get() ignored its default,
    setdefault() returned None and stored nothing, and pop() surfaced the C++
    "key does not exist" error instead of KeyError or the default.
    """
    from pyalps import ngs

    p = ngs.params({"a": 1})

    assert p["absent"] is None          # the preserved legacy quirk
    assert p.get("a") == 1
    assert p.get("absent") is None
    assert p.get("absent", 9) == 9

    assert p.setdefault("a", 5) == 1 and p["a"] == 1
    assert p.setdefault("new", 4) == 4
    assert "new" in p and p["new"] == 4

    assert p.pop("new") == 4 and "new" not in p
    assert p.pop("absent", 7) == 7
    with pytest.raises(KeyError):
        p.pop("absent")

    # observables and results raise KeyError natively, so get() is fine there,
    # but pop() needs the same replacement: MutableMapping.pop reads
    # self._MutableMapping__marker, which a copied method cannot resolve.
    obs = ngs.observables()
    obs.createRealObservable("x")
    assert obs.get("absent", 3) == 3
    assert obs.pop("absent", 3) == 3
    with pytest.raises(KeyError):
        obs.pop("absent")
    assert obs.pop("x") is not None and "x" not in obs


def test_archive_errors_use_the_typed_hierarchy(tmp_path):
    """pyalps.hdf5's exception classes must actually be raised.

    nanobind compiles extensions with -fvisibility=hidden, so the catch
    clauses in the exception translator could not match the exceptions libalps
    threw: every archive failure arrived as a bare RuntimeError carrying the
    whole ALPS_STACKTRACE, and ArchiveNotFound/ArchiveClosed were never seen.
    Assert on the message length too -- a translated exception is trimmed to
    its first line, so a multi-line message means the translator was bypassed.
    """
    import pyalps.hdf5 as hdf5

    with pytest.raises(hdf5.ArchiveNotFound) as missing:
        hdf5.archive(str(tmp_path / "does-not-exist.h5"), "r")
    assert len(str(missing.value).splitlines()) == 1

    archive = hdf5.archive(str(tmp_path / "a.h5"), "w")
    archive["/v"] = 1
    archive.close()
    with pytest.raises(hdf5.ArchiveClosed) as closed:
        archive["/v"]
    assert len(str(closed.value).splitlines()) == 1

    # every one of them derives from ArchiveError, so callers can catch broadly
    for cls in (hdf5.ArchiveNotFound, hdf5.ArchiveClosed, hdf5.InvalidPath,
                hdf5.PathNotFound, hdf5.WrongType):
        assert issubclass(cls, hdf5.ArchiveError)


def test_complex_params_hdf5_roundtrip(tmp_path):
    """Complex parameters must survive a checkpoint.

    Two separate defects made this fail. archive::set_complex() did not
    resolve its path against the current context, so the marker attribute for
    a value written at the empty path landed on the root group; and
    paramvalue::load() sent complex scalars into the vector branch, because a
    complex scalar has is_scalar() == false (it is stored as a trailing
    dimension of two reals). Rank distinguishes them: 1 for a scalar, 2 for a
    vector of any length.
    """
    from pyalps import ngs

    cases = {"scalar": 1 + 2j, "vector": [1 + 2j, 3 + 4j], "one": [5 + 6j]}
    for name, value in cases.items():
        path = str(tmp_path / ("complex-%s.h5" % name))
        with hdf5_archive(path, "w") as archive:
            ngs.params({name: value}).save(archive)
        loaded = ngs.params()
        with hdf5_archive(path, "r") as archive:
            loaded.load(archive, "/")
        got = list(loaded[name]) if isinstance(value, list) else loaded[name]
        assert got == value, "%s: %r != %r" % (name, got, value)


def test_mcbase_base_save_is_not_virtual(tmp_path):
    """Calling the base save() from an override must not re-enter the override.

    save/load were bound as pointers-to-member, which dispatch through the
    vtable, so ngs.mcbase.save(self, ar) -- and super().save(ar) -- landed back
    in the Python override and ran its body twice.
    """
    from pyalps import ngs

    class Base(ngs.mcbase):
        def __init__(self, parms):
            ngs.mcbase.__init__(self, parms, 42)
            self.measurements.createRealObservable("E")
            self.steps = 0

        def update(self):
            self.steps += 1

        def measure(self):
            self.measurements["E"] << 1.0

        def fraction_completed(self):
            return self.steps / 5.0

    for label, use_super in (("explicit", False), ("super", True)):
        calls = []

        class Override(Base):
            def save(self, archive):
                calls.append(label)
                if use_super:
                    super().save(archive)
                else:
                    ngs.mcbase.save(self, archive)

        simulation = Override({"SWEEPS": 5, "THERMALIZATION": 0, "SEED": 1})
        simulation.run(lambda: False)
        with hdf5_archive(str(tmp_path / ("mcbase-%s.h5" % label)), "w") as archive:
            simulation.save(archive)
        assert calls == [label], "%s: save() ran %d times" % (label, len(calls))


def test_params_native_bool_vector_hdf5_roundtrip(tmp_path):
    from pyalps import hdf5, ngs

    filename = str(tmp_path / "bool-params.h5")
    original = ngs.params({"flags": [True, False, True]})
    with hdf5.archive(filename, "w") as archive:
        original.save(archive)
    with hdf5.archive(filename, "r") as archive:
        loaded = ngs.params(archive, "/")
    assert loaded["flags"] == [True, False, True]


def test_observable_lshift_chains():
    from pyalps import ngs

    observables = ngs.observables()
    observables.createRealObservable("chain")
    observable = observables["chain"]
    returned = (observable << 1.0) << 2.0
    assert returned is observable
    assert ngs.observable2result(observable).count == 2


def test_standalone_observables_accept_samples():
    """ngs.createRealObservable() handles must accept measurements.

    They stopped doing so under nanobind: extensions are compiled with
    -fvisibility=hidden, so instantiating a libalps class template inside a
    binding TU emits a hidden vtable/type_info that cannot merge with
    libalps' copy, and the dynamic_cast<RecordableObservable<T>*> in
    Observable::add then fails with "Cannot add measurement to observable".
    The fix keeps construction on the libalps side; this pins it.
    """
    import numpy as np

    from pyalps import ngs

    scalar = ngs.createRealObservable("Energy")
    scalar << 1.0
    scalar << 2.0

    vector = ngs.createRealVectorObservable("Correlations")
    vector << np.array([1.0, 2.0, 3.0])

    # the container-held equivalents must keep working too
    observables = ngs.observables()
    observables.createRealObservable("Energy")
    observables["Energy"] << 1.5


def test_observables_item_deletion():
    from pyalps import ngs

    observables = ngs.observables()
    observables.createRealObservable("a")
    observables.createRealObservable("b")
    del observables["a"]
    assert "a" not in observables and "b" in observables
    observables.clear()
    assert len(observables) == 0


def test_mapping_views_are_set_like():
    """keys/values/items must be MutableMapping views, not one-shot iterators.

    Boost.Python's map_indexing_suite defined none of the three, so on the
    legacy build they resolved through MutableMapping to KeysView/ValuesView/
    ItemsView: sized, re-iterable and set-like. The nanobind port must keep
    that, which means NOT defining them natively in C++ -- pyalps/ngs.py only
    grafts a mixin onto names the extension type leaves alone.
    """
    from collections.abc import MutableMapping

    from pyalps import ngs

    observables = ngs.observables()
    observables.createRealObservable("a")
    observables.createRealObservable("b")

    for mapping in (observables, ngs.params({"a": 1, "b": 2})):
        keys = mapping.keys()
        # sized, and re-iterable (a nanobind iterator is exhausted after one pass)
        assert len(keys) == 2
        assert sorted(keys) == ["a", "b"]
        assert sorted(keys) == ["a", "b"]
        # set-like
        assert keys & {"a"} == {"a"}
        assert keys | {"c"} == {"a", "b", "c"}

        items = mapping.items()
        assert len(items) == 2
        assert sorted(k for k, _ in items) == ["a", "b"]
        assert sorted(k for k, _ in items) == ["a", "b"]

        values = mapping.values()
        assert len(values) == 2
        assert len(list(values)) == 2
        assert len(list(values)) == 2

    # `results` is the third mapping type and goes through the same shim, but
    # it is deliberately not constructible from Python -- master bound it with
    # boost::python::no_init and the port binds no nb::init<> either -- so the
    # view semantics are asserted here only through the two types that are.
    for _name in ("keys", "values", "items"):
        assert getattr(ngs.results, _name) is getattr(MutableMapping, _name), (
            "results.%s must come from the MutableMapping mixin, not a native "
            "one-shot nanobind iterator" % _name
        )


def test_mcbase_save_load_overrides_reach_cpp_dispatch():
    from pyalps import ngs
    from pyalps.cxx import pyngshdf5_c

    calls = []

    class Simulation(ngs.mcbase):
        def update(self):
            pass

        def measure(self):
            pass

        def fraction_completed(self):
            return 1.0

        def save(self, archive):
            calls.append("save")
            super().save(archive)

        def load(self, archive):
            calls.append("load")
            super().load(archive)

    simulation = Simulation({"SEED": 42})
    # the base save/load expects a non-empty measurements container
    simulation.measurements << ngs.RealObservable("energy")
    simulation.measurements["energy"] << 1.0
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "checkpoint.h5")
        # Drive a real C++-side checkpoint rather than calling the base
        # binding: `archive[path] = simulation` hands the object to the C++
        # save path, which must reach the Python override. Calling
        # ngs.mcbase.save(simulation, archive) would NOT test this -- that is
        # the base implementation and deliberately does not dispatch
        # virtually, so it cannot re-enter the override (see
        # test_mcbase_base_save_is_not_virtual).
        archive = pyngshdf5_c.hdf5_archive_impl(path, "w")
        archive["/simulation"] = simulation
        del archive
        assert calls == ["save"]
        # the override's super().save() must have written the real payload
        archive = pyngshdf5_c.hdf5_archive_impl(path, "r")
        assert "measurements" in archive.list_children("/simulation")
        archive.set_context("/simulation")
        simulation.load(archive)
        del archive
        assert calls == ["save", "load"]


def test_accumulator_result_inplace_identity():
    from pyalps.cxx.pyngsaccumulator_c import error_accumulator

    accumulator = error_accumulator()
    accumulator(1.0)
    accumulator(2.0)
    result = accumulator.result()
    alias = result
    alias += 1.0
    assert alias is result
    assert np.isclose(result.mean(), 2.5)


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
        test_params_from_parameter_file,
        test_alea_numpy_and_mcdata_operators,
        test_alea_mcanalyze_surface,
        test_packaged_xml_stylesheets,
        test_ngs_observable_containers,
        test_name_encoding_roundtrip,
        test_accumulator_surface,
        test_optional_application_extension_surface,
        test_params_mapping_equality_and_value_ladder,
        test_params_mapping_mixins_handle_none_getitem,
        test_observable_lshift_chains,
        test_standalone_observables_accept_samples,
        test_observables_item_deletion,
        test_mapping_views_are_set_like,
        test_mcbase_save_load_overrides_reach_cpp_dispatch,
        test_accumulator_result_inplace_identity,
    ):
        test()
    print("pyalps binding surface: green")
