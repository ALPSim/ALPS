# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/large_active_space.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 1fe88978ad459480ed34dab9a7ac9738bdb8a387b8c80a2c0f85538b8de9b5d0.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""Get a usable wavefunction past the CI-reconstruction limit.

Requires a QDK build with MPSContainer (PR #751). At 14 orbitals half filled
there are 11.8 million determinants: ~94 MB of float64 amplitudes, but
~12 GB at the classic adapter's conservative 1 kB per determinant estimate. The calculator
returns an MPS and the active-space RDMs instead. See rdm_observables.py for the
ALPS-native demonstration that also runs with released QDK packages.
"""

import numpy as np
from qdk_chemistry.algorithms import create
from qdk_chemistry.data import LatticeGraph
from qdk_chemistry.utils.model_hamiltonians import create_hubbard_hamiltonian

import qdk_chemistry_alps  # noqa: F401 — register ALPS calculators

n_orbitals = 14
hamiltonian = create_hubbard_hamiltonian(
    LatticeGraph.chain(n_orbitals), epsilon=0.0, t=1.0, U=4.0
)

calc = create("multi_configuration_calculator", "alps_classic_dmrg")
calc.settings().set("max_bond_dim", 64)
calc.settings().set("wavefunction_container", "mps")
calc.settings().set("sweeps", 8)
energy, wavefunction = calc.run(hamiltonian, 7, 7)

print(f"container : {wavefunction.get_container_type()}")
print(f"energy    : {energy:.6f}")

container = wavefunction.get_container()
print("method    : alps_classic_dmrg")
print(f"bond dim  : {container.max_bond_dimension}")

# The stored RDMs reproduce the energy through the integrals.
dm1 = np.asarray(wavefunction.get_active_one_rdm_spin_traced())
dm2 = np.asarray(wavefunction.get_active_two_rdm_spin_traced())
dm2 = dm2.reshape(n_orbitals, n_orbitals, n_orbitals, n_orbitals)

h1 = np.asarray(hamiltonian.get_one_body_integrals()[0])
two = hamiltonian.get_two_body_integrals()
eri = np.asarray(two[0] if isinstance(two, (tuple, list)) else two)
eri = eri.reshape(n_orbitals, n_orbitals, n_orbitals, n_orbitals)

from_rdms = (
    np.einsum("pq,pq->", h1, dm1)
    + 0.5 * np.einsum("pqrs,pqrs->", eri, dm2)
    + float(hamiltonian.get_core_energy())
)
print(f"energy from stored RDMs : {abs(from_rdms - energy):.2e}  (difference)")

# Natural occupations are eigenvalues of the stored spin-summed 1-RDM.
# This also works with QDK builds whose MPS occupation accessor is absent.
natural = np.linalg.eigvalsh(dm1)[::-1]
print(f"natural occupations sum : {natural.sum():.6f}  (= electron count)")

assert wavefunction.get_container_type() == "mps"
assert abs(from_rdms - energy) < 1e-7
assert abs(natural.sum() - 14) < 1e-7
