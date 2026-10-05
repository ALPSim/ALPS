#!/usr/bin/env python3
"""Run the nanobind-exported C++ simulation from Python."""

import pyalps.ngs as ngs

import ising_c


simulation = ising_c.sim(ngs.params({"SEED": 7, "SWEEPS": 10}))
simulation.run(lambda: False)
results = simulation.collectResults()
print(results)
