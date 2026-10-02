#!/usr/bin/env python3
"""Exercise a simulation written entirely in Python against pyalps.ngs.

Unlike 7_python_extend, `ising.sim` here does not derive from ngs.mcbase: it
composes ngs.params, ngs.random01 and ngs.createRealObservable by hand and
implements the scheduler protocol (update / measure / fraction_completed / run)
in Python. Its save() and load() are the archive-object protocol -- pyalps
finds a Python-defined save() on the object and calls it with the archive.

This tutorial had no test of any kind, and both halves of its checkpoint were
broken: save() could not store an ngs.params or the measurements dict through
`archive[path] = ...`, and load() wrote to the archive instead of reading from
it, passed a value where params.load wants an archive, and called double().
"""

import os
import tempfile

import numpy as np

import ising
import pyalps.hdf5 as hdf5


PARAMETERS = {"L": 12, "THERMALIZATION": 5, "SWEEPS": 20, "T": 2.0, "SEED": 42}

OBSERVABLES = ["Correlations", "Energy", "Magnetization",
               "Magnetization^2", "Magnetization^4"]


simulation = ising.sim(PARAMETERS)

assert int(simulation.parameters["L"]) == 12
assert sorted(simulation.result_names()) == OBSERVABLES
assert simulation.fraction_completed() == 0.0

# run() returns True when it finished rather than being stopped.
assert simulation.run(lambda: False)
assert simulation.fraction_completed() >= 1.0

results = simulation.collectResults()
assert sorted(results) == OBSERVABLES
for name in OBSERVABLES:
    assert results[name].count == PARAMETERS["SWEEPS"], name
assert -1.0 <= results["Energy"].mean <= 1.0
assert len(results["Correlations"].mean) == PARAMETERS["L"]

with tempfile.TemporaryDirectory() as directory:
    checkpoint = os.path.join(directory, "ising.clone0.h5")

    # `archive['/'] = simulation` dispatches to the Python-defined save(),
    # which in turn stores an ngs.params and a dict of observables.
    with hdf5.archive(checkpoint, "w") as archive:
        archive["/"] = simulation

    with hdf5.archive(checkpoint, "r") as archive:
        assert sorted(archive.list_children("/parameters")) == [
            "L", "SEED", "SWEEPS", "T", "THERMALIZATION"]
        clone = "/simulation/realizations/0/clones/0"
        assert sorted(archive.list_children(clone + "/measurements")) == OBSERVABLES
        assert sorted(archive.list_children(clone + "/checkpoint")) == [
            "engine", "spins", "sweeps"]

    restored = ising.sim(PARAMETERS)
    with hdf5.archive(checkpoint, "r") as archive:
        restored.load(archive)

    assert restored.sweeps == simulation.sweeps
    np.testing.assert_array_equal(restored.spins, simulation.spins)
    assert int(restored.parameters["L"]) == int(simulation.parameters["L"])

    # The restored measurements carry the samples from before the checkpoint.
    restored_results = restored.collectResults()
    for name in OBSERVABLES:
        assert restored_results[name].count == results[name].count, name

    # A reloaded RNG must continue the original stream, not restart it.
    assert [restored.random() for _ in range(4)] == \
           [simulation.random() for _ in range(4)]

print("native python simulation: ok")
