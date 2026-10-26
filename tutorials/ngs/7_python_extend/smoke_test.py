#!/usr/bin/env python3
"""Exercise a Python subclass of ngs.mcbase driving the C++ scheduler.

This is the migration's headline capability: `ising.sim` derives from
ngs.mcbase, so C++ calls back into Python for update / measure /
fraction_completed through nanobind's trampoline, while parameters, the RNG and
the measurement container come from the C++ base.

Its save()/load() also exercise cooperative checkpointing -- the subclass calls
ngs.mcbase.save(self, ar) for the base state and then adds its own. That only
works because the base save is bound with a qualified call: dispatching it
through the vtable would re-enter the subclass override and run the body twice.

The tutorial had no test, and load() had never run: it called double() and read
its two values from paths that save() does not write.
"""

import os
import tempfile

import numpy as np

import ising
import pyalps.hdf5 as hdf5
import pyalps.ngs as ngs


PARAMETERS = {"L": 12, "THERMALIZATION": 5, "SWEEPS": 20, "T": 2.0}

OBSERVABLES = ["Correlations", "Energy", "Magnetization",
               "Magnetization^2", "Magnetization^4"]

CLONE = "/simulation/realizations/0/clones/0"


simulation = ising.sim(PARAMETERS, 42)

# The subclass really is an mcbase, and the base descriptors serve it.
assert isinstance(simulation, ngs.mcbase)
assert int(simulation.parameters["L"]) == 12
assert sorted(simulation.measurements) == OBSERVABLES
assert 0.0 <= simulation.random() < 1.0
assert simulation.fraction_completed() == 0.0

# run() drives update()/measure() in Python from the C++ scheduler.
assert simulation.run(lambda: False)
assert simulation.fraction_completed() >= 1.0
assert simulation.sweeps == PARAMETERS["THERMALIZATION"] + PARAMETERS["SWEEPS"]

results = ngs.collectResults(simulation)
assert sorted(results) == OBSERVABLES
for name in OBSERVABLES:
    assert results[name].count == PARAMETERS["SWEEPS"], name
assert -1.0 <= results["Energy"].mean <= 1.0
assert len(results["Correlations"].mean) == PARAMETERS["L"]

# The energy must match the exact periodic chain, not just lie in [-1, 1]: a
# sampler with an inverted Metropolis acceptance gives the same magnitude with
# the opposite sign. At T = 1 and L = 12, <E>/L = -(t + t^(L-1)) / (1 + t^L)
# with t = tanh(1), i.e. -0.782; 1000 sweeps keep the statistical error well
# inside the tolerance.
PHYSICS_PARAMETERS = {"L": 12, "THERMALIZATION": 100, "SWEEPS": 1000, "T": 1.0}
coupling = np.tanh(1.0 / PHYSICS_PARAMETERS["T"])
length = PHYSICS_PARAMETERS["L"]
exact_energy = -(coupling + coupling**(length - 1)) / (1 + coupling**length)

physics = ising.sim(PHYSICS_PARAMETERS, 42)
assert physics.run(lambda: False)
energy = ngs.collectResults(physics)["Energy"].mean
assert abs(energy - exact_energy) < 0.2, (energy, exact_energy)

# A stop callback that fires immediately must report "stopped", not "finished".
stopped = ising.sim(PARAMETERS, 42)
assert not stopped.run(lambda: True)

with tempfile.TemporaryDirectory() as directory:
    checkpoint = os.path.join(directory, "ising.clone0.h5")

    with hdf5.archive(checkpoint, "w") as archive:
        archive["/"] = simulation

    with hdf5.archive(checkpoint, "r") as archive:
        # ngs.mcbase.save writes relative to the archive's context, which is
        # "/" here, so the base state lands at the root -- while the subclass
        # writes its own two values at absolute paths under the clone. The
        # asymmetry is the tutorial's; load() reads each back from where it
        # was written. mcbase defines SEED, with its default of 42.
        assert sorted(archive.list_children("/parameters")) == [
            "L", "SEED", "SWEEPS", "T", "THERMALIZATION"]
        assert sorted(archive.list_children("/measurements")) == OBSERVABLES
        # The random number generator's state, saved as one string.
        assert archive.is_data("/checkpoint/engine")
        assert int(archive[CLONE + "/checkpoint/sweeps"]) == simulation.sweeps
        np.testing.assert_array_equal(
            archive[CLONE + "/checkpoint/spins"], simulation.spins)

    restored = ising.sim(PARAMETERS, 42)
    with hdf5.archive(checkpoint, "r") as archive:
        restored.load(archive)

    assert restored.sweeps == simulation.sweeps
    np.testing.assert_array_equal(restored.spins, simulation.spins)
    assert restored.length == simulation.length
    assert restored.beta == simulation.beta

    restored_results = ngs.collectResults(restored)
    assert sorted(restored_results) == OBSERVABLES
    for name in OBSERVABLES:
        assert restored_results[name].count == results[name].count, name

    # Results also go out through the public writer the tutorial uses.
    measurements = os.path.join(directory, "ising.h5")
    with hdf5.archive(measurements, "w") as archive:
        ngs.saveResults(restored_results, restored.parameters,
                        archive, "/simulation/results")
    with hdf5.archive(measurements, "r") as archive:
        assert sorted(archive.list_children("/simulation/results")) == OBSERVABLES

print("python subclass of mcbase: ok")
