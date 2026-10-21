# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/converging_dmrg.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 72319cd57176802de1ff4b5c63c400be5b8bbf1cf1c47b4aa799db0905b28567.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""Check that a DMRG run is converged, rather than assuming it.

This scan checks energy stability under increasing bond dimension. Also
check sweeps, Hamiltonian tolerances and the observables of interest. It
uses the optional MPSContainer build; omit the mps setting to run this small
problem through the released QDK environment's default CI-return path.
"""

from qdk_chemistry.algorithms import create
from qdk_chemistry.data import LatticeGraph
from qdk_chemistry.utils.model_hamiltonians import create_hubbard_hamiltonian

import qdk_chemistry_alps  # noqa: F401 — register ALPS calculators

hamiltonian = create_hubbard_hamiltonian(
    LatticeGraph.chain(8), epsilon=0.0, t=1.0, U=8.0
)

energies = {}
for bond_dimension in [4, 8, 16, 32, 64, 128, 256]:
    calc = create("multi_configuration_calculator", "alps_classic_dmrg")
    calc.settings().set("max_bond_dim", bond_dimension)
    calc.settings().set("wavefunction_container", "mps")
    energies[bond_dimension], _ = calc.run(hamiltonian, 4, 4)

best = min(energies.values())
print(f"{'bond dim':>9}  {'energy':>14}  {'above best':>12}")
for bond_dimension, energy in energies.items():
    print(f"{bond_dimension:>9}  {energy:>14.9f}  {energy - best:>12.2e}")

# More states need not monotonically improve a finite-sweep optimization.
# Show observed differences and compare with an independent small-sector solve.
reference, _ = create("multi_configuration_calculator", "macis_cas").run(hamiltonian, 4, 4)
print(f"MACIS reference: {reference:.12f}; final-state error: {energies[256]-reference:.3e}")
assert abs(energies[256]-reference) < 1e-7
