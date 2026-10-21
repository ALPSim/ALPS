# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/quickstart.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 c557f150d0ad9e8f753c6326122883bb8c401ebca22f0f8b64a7a0cb5465885e.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""Solve a Hubbard chain with both backends and compare them.

The shortest end-to-end path: build a model with qdk-chemistry's own
factory, solve it two ways, and check they agree.
"""

from qdk_chemistry.algorithms import create
from qdk_chemistry.data import LatticeGraph
from qdk_chemistry.utils.model_hamiltonians import create_hubbard_hamiltonian

import qdk_chemistry_alps  # noqa: F401 — register ALPS calculators

# A six-site Hubbard chain at half filling.
hamiltonian = create_hubbard_hamiltonian(
    LatticeGraph.chain(6), epsilon=0.0, t=1.0, U=4.0
)

# Exact diagonalization in the full fixed-particle-number sector.
ed = create("multi_configuration_calculator", "alps_ed")
energy_ed, wavefunction = ed.run(hamiltonian, 3, 3)

# DMRG: optimize an MPS, then check its accuracy against ED here.
dmrg = create("multi_configuration_calculator", "alps_classic_dmrg")
dmrg.settings().set("max_bond_dim", 128)
energy_dmrg, _ = dmrg.run(hamiltonian, 3, 3)

print(f"exact diagonalization : {energy_ed:.10f}")
print(f"DMRG                  : {energy_dmrg:.10f}")
print(f"difference            : {abs(energy_ed - energy_dmrg):.2e}")
print(f"determinants returned : {len(wavefunction.get_coefficients())}")

assert abs(energy_ed - energy_dmrg) < 1e-8
