# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/custom_spin_export.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 4221ed3609eed3aaf261390a93c95b3825e4034b02ccc85555ac3cd54f5813b0.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""ALPS → QDK: a cluster–Ising model with a three-site term and a Y field."""

import time

from pyalps import models
import numpy as np
from demo_support import model_matrix
from pyalps.models import Hamiltonian, spin as Spin
from demo_support import bit_reversal, finish, parser, pauli_matrix, site_matrix

from qdk_chemistry_alps.interop import hamiltonian_to_qubit_operator


def model(n=6):
    spin = Spin()
    h = Hamiltonian(models.lattice.chain(n, bc="open"), spin)
    reference = .13 * np.eye(2**n, dtype=complex)
    h.add(.13 * models.located(spin, "I", 0))
    for i in range(1, n - 1):
        h.add(-8 * models.located(spin, "Sz", i - 1) * models.located(spin, "Sx", i) * models.located(spin, "Sz", i + 1))
        reference -= site_matrix(n, {i - 1: "Z", i: "X", i + 1: "Z"})
    fields = np.linspace(-.19, .27, n)
    for i, field in enumerate(fields):
        h.add(-.74 * models.located(spin, "Sx", i) + 2 * field * models.located(spin, "Sz", i))
        reference += -.37 * site_matrix(n, {i: "X"}) + field * site_matrix(n, {i: "Z"})
    h.add(.22 * models.located(spin, "Sy", 1))
    reference += .11 * site_matrix(n, {1: "Y"})
    return h, reference


def plot(plt, data):
    figure, axes = plt.subplots(1, 2, figsize=(9, 3.5), layout="constrained")
    axes[0].plot(data["reference_spectrum"], ".", color="#087e8b", label="Independent reference")
    axes[0].plot(data["exported_spectrum"], "+", color="#d68c22", label="QDK Pauli operator")
    axes[0].set(xlabel="Eigenvalue index", ylabel="Energy / J", title="Hamiltonian spectra")
    axes[0].legend()
    axes[1].bar(["Pauli Z₀", "Pauli Y₁"], data["reference_expectations"], color="#087e8b", label="ALPS product state")
    axes[1].plot([0, 1], data["qubit_expectations"], "x", color="#c44e52", markersize=10, label="Mapped QDK state")
    axes[1].set(ylabel="Expectation", title="Convention-sensitive observables")
    axes[1].legend()
    return figure


def main():
    args = parser(__doc__).parse_args()
    started = time.perf_counter()
    h, exact = model()
    operator = hamiltonian_to_qubit_operator(h)
    exported = pauli_matrix(operator)
    reverse = bit_reversal(6)
    np.testing.assert_allclose(model_matrix(h), exact, atol=1e-12, rtol=0)
    np.testing.assert_allclose(exported[np.ix_(reverse, reverse)], exact, atol=1e-12, rtol=0)
    np.testing.assert_allclose(exported, exported.conj().T, atol=1e-12, rtol=0)
    assert np.max(abs(exported - exact)) > .1  # wrong-endianness negative control
    # Distinct complex local states expose both reversal and the sign of Y.
    psi = np.ones(1, dtype=complex)
    for i in range(6):
        theta, phase = .17 + .09 * i, .23 + .21 * i
        psi = np.kron(psi, [np.cos(theta), np.exp(1j * phase) * np.sin(theta)])
    mapped = psi[reverse]
    spin = h.site_for(0)
    expected_values, actual_values = [], []
    for expression, factors in [(2 * models.located(spin, "Sz", 0), {0: "Z"}), (2 * models.located(spin, "Sy", 1), {1: "Y"})]:
        observable = Hamiltonian(h.lattice, spin).add(expression)
        expected_values.append(np.vdot(psi, site_matrix(6, factors) @ psi).real)
        actual_values.append(np.vdot(mapped, pauli_matrix(hamiltonian_to_qubit_operator(observable)) @ mapped).real)
    np.testing.assert_allclose(actual_values, expected_values, atol=1e-12, rtol=0)
    np.testing.assert_allclose(exported @ mapped, (exact @ psi)[reverse], atol=1e-12, rtol=0)
    words = operator.to_dict()["pauli_strings"]
    print(f"Cluster–Ising: 6 qubits, {len(words)} Pauli terms, maximum support {max(6 - w.count('I') for w in words)}")
    print("Full matrix, constant, complex Y term, state action and mapped observables agree.")
    print("Static model export; this example does not evolve the three-site interaction.")
    finish(args, {"model": {"sites": 6, "J": 1., "h": .37, "longitudinal_fields": np.linspace(-.19, .27, 6),
                            "y_field_site": 1, "y_field": .11, "constant": .13, "boundary": "open"},
                  "pauli_operator": operator.to_dict(), "matrix_error": np.max(abs(exported[np.ix_(reverse, reverse)] - exact)),
                  "reference_spectrum": np.linalg.eigvalsh(exact), "exported_spectrum": np.linalg.eigvalsh(exported),
                  "reference_expectations": expected_values, "qubit_expectations": actual_values,
                  "basis_map": "site 0 most significant in ALPS; qubit 0 least significant in QDK Pauli matrix"}, started, plot)


if __name__ == "__main__":
    main()
