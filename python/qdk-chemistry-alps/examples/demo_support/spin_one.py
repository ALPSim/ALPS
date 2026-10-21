# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/demo_support/spin_one.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 68749a20564fffa0525291c860ab3011d6185c220982fdddb3e8b2246a54da8d.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""Four-site spin-1 models and small independent matrix references."""

import numpy as np
from pyalps import models
from pyalps.models.lattice import chain
from pyalps.models import Hamiltonian, spin as Spin

COUPLINGS = np.array([1., .8, 1.2])
FIELDS = np.array([.07, -.11, .16, -.03])


def embed(matrix, site, n=4):
    result = np.ones((1, 1), dtype=complex)
    for i in range(n):
        result = np.kron(result, matrix if i == site else np.eye(3))
    return result


def model(anisotropy, *, biquadratic=.18, transverse=0., y_field=0.):
    spin = Spin(1.)
    h = Hamiltonian(chain(4, bc="open"), spin,
                    conserved_quantum_numbers=["Sz"] if transverse == y_field == 0 else [])
    # Explicit angular momentum matrices, independent of ALPS operators.
    x = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]) / np.sqrt(2)
    y = np.array([[0, -1j, 0], [1j, 0, -1j], [0, 1j, 0]]) / np.sqrt(2)
    z = np.diag([1., 0., -1.])
    axes = [[embed(s, i) for i in range(4)] for s in (x, y, z)]
    reference = .37 * np.eye(81, dtype=complex)
    h.add(.37 * models.located(spin, "I", 0))
    for i, coupling in enumerate(COUPLINGS):
        expression = (models.located(spin, "Sz", i)*models.located(spin, "Sz", i+1) + .5*models.located(spin, "Sp", i)*models.located(spin, "Sm", i+1)
                      + .5*models.located(spin, "Sm", i)*models.located(spin, "Sp", i+1))
        h.add(coupling*expression + biquadratic*expression*expression)
        exchange = sum(axis[i] @ axis[i+1] for axis in axes)
        reference += coupling*exchange + biquadratic*exchange@exchange
    for i, field in enumerate(FIELDS):
        h.add(anisotropy*models.located(spin, "Sz", i)*models.located(spin, "Sz", i) + field*models.located(spin, "Sz", i))
        reference += anisotropy*axes[2][i]@axes[2][i] + field*axes[2][i]
    if transverse:
        h.add(transverse*models.located(spin, "Sx", 1))
        reference += transverse*axes[0][1]
    if y_field:
        h.add(y_field*models.located(spin, "Sy", 2))
        reference += y_field*axes[1][2]
    return h, reference, axes


def reference_indices(scheme):
    """Independent declared mapping: physical site 0 is most significant."""
    indices = []
    for local in np.ndindex(3, 3, 3, 3):
        codes = [j if scheme == "binary" else j ^ (j >> 1) for j in local]
        indices.append(sum(code << (2*i) for i, code in enumerate(codes)))
    return np.asarray(indices)


def qubit_matrix(operator):
    """Read QDK terms and apply Pauli actions to register basis vectors."""
    data = operator.to_dict()
    coefficients = np.asarray(data["coefficients"]["real"]) + 1j*np.asarray(data["coefficients"]["imag"])
    columns = np.arange(2**operator.num_qubits)
    result = np.zeros((len(columns), len(columns)), dtype=complex)
    for word, coefficient in zip(data["pauli_strings"], coefficients, strict=True):
        rows = columns.copy()
        values = np.full(len(columns), coefficient, dtype=complex)
        for q, letter in enumerate(word[::-1]):
            bits = (columns >> q) & 1
            if letter in "XY":
                rows ^= 1 << q
            if letter == "Z":
                values *= 1 - 2*bits
            if letter == "Y":
                values *= 1j*(1 - 2*bits)
        result[rows, columns] += values
    return result
