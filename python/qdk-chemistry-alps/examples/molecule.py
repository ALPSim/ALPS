# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/molecule.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 fcfcf356d4bf3ae8c24abb7a791b36c2ce8c019bc582280e12a5580faf5edc5b.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""Run a molecular active space from PySCF integrals.

Needs pyscf. Other sources can supply real restricted integrals (h1, eri)
in chemist notation, expressed in an orthonormal one-particle basis.
"""

from pyscf import ao2mo, fci, gto, scf

from qdk_chemistry_alps.abinitio import build_abinitio_hamiltonian, solve_ground_state

mol = gto.M(atom="Li 0 0 0; H 0 0 1.6", basis="sto-3g", verbose=0)

mean_field = scf.RHF(mol)
mean_field.conv_tol = 1e-12
mean_field.kernel()

# Transform to the molecular-orbital basis.
coefficients = mean_field.mo_coeff
n_orbitals = coefficients.shape[1]
h1 = coefficients.T @ mean_field.get_hcore() @ coefficients
eri = ao2mo.restore(1, ao2mo.kernel(mol, coefficients), n_orbitals)
e_core = float(mol.energy_nuc())

# Build the ALPS operator and solve it in the (N, Sz) sector.
H = build_abinitio_hamiltonian(h1, eri)
result = solve_ground_state(H, n_alpha=2, n_beta=2)
energy = float(result.eigenvalues[0]) + e_core

reference = fci.direct_spin1.FCI().kernel(h1, eri, n_orbitals, (2, 2))[0] + e_core

print(f"orbitals        : {n_orbitals}")
print(f"Hartree-Fock    : {mean_field.e_tot:.10f}")
print(f"ALPS ED         : {energy:.10f}")
print(f"PySCF FCI       : {reference:.10f}")
print(f"difference      : {abs(energy - reference):.2e}")
print(f"correlation     : {energy - mean_field.e_tot:.6f}")

assert abs(energy - reference) < 1e-8
