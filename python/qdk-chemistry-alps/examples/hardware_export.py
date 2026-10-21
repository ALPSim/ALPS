# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/hardware_export.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 0410474084abc9cff25961fbe43e2124f170056bf446e346c5658fe0989834fc.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""Export one fermionic Hamiltonian as a QDK qubit operator.

The ALPS and QDK representations must describe the same quantum model.
This example checks their full spectra by numerical diagonalization on an
ordinary computer; it does not construct a circuit or run on a device.
Spectral agreement alone does not establish a correct state/basis map.
"""

from pyalps import models
import numpy as np
from demo_support import model_matrix
from pyalps.models import hubbard

from qdk_chemistry_alps import interop

# A fermionic model built on the ALPS side.
H = hubbard(models.lattice.chain(3, bc="open"), t=1.0, U=4.0, mu=0.7)

for encoding in ["jordan-wigner", "bravyi-kitaev", "parity"]:
    qubit_op = interop.fermionic_hamiltonian_to_qubit_operator(H, encoding)

    # Compare the full spectra, not just the ground state: the operators
    # must agree on every state, constant term included.
    fock = np.sort(np.linalg.eigvalsh(model_matrix(H)))
    as_qubits = interop.qubit_operator_to_hamiltonian(qubit_op)
    qubit = np.sort(np.linalg.eigvalsh(model_matrix(as_qubits)))

    terms = len(qubit_op.to_dict()["pauli_strings"])
    print(
        f"{encoding:18s} {qubit_op.num_qubits} qubits, {terms:2d} Pauli terms, "
        f"spectrum agrees to {np.max(np.abs(fock - qubit)):.1e}"
    )
    assert np.max(np.abs(fock - qubit)) < 1e-12

# Metadata travels with the operator, so downstream code knows the layout.
meta = interop.fermionic_hamiltonian_to_qubit_operator(H).to_dict()
print(f"encoding={meta['encoding']!r} mode_order={meta['fermion_mode_order']!r}")
