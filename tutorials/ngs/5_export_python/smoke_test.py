#!/usr/bin/env python3
"""Exercise the public downstream simulation-export compatibility helper."""

import os
import tempfile

# Importing the consumer first verifies its wheel-runtime rpath. Its module
# initializer loads the owning pyalps bindings before registering C++ types.
import ising_c
import pyalps.hdf5 as hdf5
import pyalps.ngs as ngs


parameters = ngs.params({"SEED": 7, "SWEEPS": 10})
simulation = ising_c.sim(parameters)

assert issubclass(ising_c.sim, ngs.mcbase)
assert isinstance(simulation, ngs.mcbase)
assert int(simulation.parameters["SWEEPS"]) == 10
assert len(simulation.measurements) == 1
assert 0.0 <= simulation.random() < 1.0
assert simulation.run(lambda: False)
assert simulation.resultNames() == ["Magnetization"]
before = simulation.collectResults()
assert before["Magnetization"].count == 10

with tempfile.TemporaryDirectory() as directory:
    checkpoint = os.path.join(directory, "ising.h5")
    with hdf5.archive(checkpoint, "w") as archive:
        simulation.save(archive)

    restored = ising_c.sim(parameters)
    with hdf5.archive(checkpoint, "r") as archive:
        restored.load(archive)

    after = restored.collectResults()
    assert restored.resultNames() == simulation.resultNames()
    assert after["Magnetization"].count == before["Magnetization"].count
    assert after["Magnetization"].mean == before["Magnetization"].mean

print("downstream nanobind export: ok")
