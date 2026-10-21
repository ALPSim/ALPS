# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/fermion_pairing.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 06becfd2b740ce60ae0024a2d0cd7ff1267b0a85b4493efe448afa3cfe036355.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""A spinless pairing chain: native ALPS → QDK without chemistry integrals."""

import time

from pyalps import models
import numpy as np
from demo_support import export_operator
from demo_support import model_matrix
from pyalps.models import Hamiltonian, spinless_fermion as SpinlessFermion
from demo_support import finish, fock_matrix, parser, pauli_matrix

from qdk_chemistry_alps.interop import pauli_operator_to_qubit_operator


def model(n=6):
    site = SpinlessFermion()
    # Pairing does NOT conserve N: deliberately leave conserved_qns empty.
    h = Hamiltonian(models.lattice.chain(n, bc="open"), site)
    products = []
    for i in range(n - 1):
        h.add(-models.located(site, "c_dag", i) * models.located(site, "c", i + 1), add_hc=True)
        h.add((.3 + .12j) * models.located(site, "c", i) * models.located(site, "c", i + 1), add_hc=True)
        h.add(.15 * models.located(site, "n", i) * models.located(site, "n", i + 1))
        products += [(-1., [(i, True), (i + 1, False)]),
                     (-1., [(i + 1, True), (i, False)]),
                     (.3 + .12j, [(i, False), (i + 1, False)]),
                     (.3 - .12j, [(i + 1, True), (i, True)]),
                     (.15, [(i, True), (i, False), (i + 1, True), (i + 1, False)])]
    for i in range(n):
        mu = .7 + .05 * i
        h.add(-mu * models.located(site, "n", i) + .5 * mu * models.located(site, "I", i))
        products += [(-mu, [(i, True), (i, False)]), (.5 * mu, [])]
    return h, fock_matrix(list(range(2**n)), products)


def plot(plt, data):
    figure, axes = plt.subplots(1, 2, figsize=(9, 3.5), layout="constrained")
    for parity, color in [("even", "#087e8b"), ("odd", "#c44e52")]:
        values = data[f"{parity}_spectrum"][:8]
        axes[0].plot(range(8), values, "o-", color=color, label=f"{parity} parity")
    axes[0].set(xlabel="Eigenvalue within parity block", ylabel="Energy / t", title="Pairing preserves parity, not particle number")
    axes[0].legend()
    axes[1].bar(range(6), data["occupations"], color="#087e8b")
    axes[1].set(xlabel="Site", ylabel="⟨nᵢ⟩", ylim=(0, 1), title="Exported-operator ground-state occupations")
    return figure


def main():
    args = parser(__doc__).parse_args()
    started = time.perf_counter()
    h, exact = model()
    encoded = export_operator(h)
    qdk = pauli_operator_to_qubit_operator(encoded)
    exported = pauli_matrix(qdk)
    indices = [encoded.state_index(row) for row in np.array(list(__import__("itertools").product(range(2), repeat=6)))]
    np.testing.assert_allclose(exported, exact, atol=1e-12, rtol=0)
    np.testing.assert_allclose(exported[np.ix_(indices, indices)], model_matrix(h), atol=1e-12, rtol=0)
    energies, vectors = np.linalg.eigh(exported)
    psi = vectors[:, 0]
    occupation = np.array([[(state >> i) & 1 for i in range(6)] for state in range(64)])
    number = occupation.sum(axis=1)
    parity = (-1.)**number
    np.testing.assert_allclose(exported * parity - parity[:, None] * exported, 0, atol=1e-12)
    assert np.max(abs(exported * number - number[:, None] * exported)) > .1
    site = h.site_for(0)
    occupations = []
    for i in range(6):
        observable = export_operator(h, observable=models.located(site, "n", i))
        occupations.append(np.vdot(psi, pauli_matrix(pauli_operator_to_qubit_operator(observable)) @ psi).real)
    np.testing.assert_allclose(occupations, abs(psi)**2 @ occupation, atol=1e-12)
    pairing = pauli_matrix(pauli_operator_to_qubit_operator(export_operator(h, observable=models.located(site, "c", 0) * models.located(site, "c", 1))))
    reference_pairing = fock_matrix(list(range(64)), [(1., [(0, False), (1, False)])])
    np.testing.assert_allclose(pairing @ psi, reference_pairing @ psi, atol=1e-12)
    pair_expectation = np.vdot(psi, pairing @ psi)
    variance = float(abs(psi)**2 @ number**2 - (abs(psi)**2 @ number)**2)
    assert variance > .01
    data = {"model": {"sites": 6, "t": 1., "pairing_real": .3, "pairing_imag": .12,
                      "mu": .7 + .05 * np.arange(6), "V": .15, "boundary": "open"},
            "layout": {"site_qubits": encoded.site_qubits, "fermion_mode_qubits": encoded.fermion_mode_qubits},
            "discarded_l1": encoded.discarded_l1, "pauli_operator": qdk.to_dict(),
            "energy": energies[0], "occupations": occupations,
            "particle_number_variance": variance,
            "pair_expectation": {"real": pair_expectation.real, "imag": pair_expectation.imag}}
    for label, sign in [("even", 1), ("odd", -1)]:
        block = np.flatnonzero(parity == sign)
        data[f"{label}_spectrum"] = np.linalg.eigvalsh(exported[np.ix_(block, block)])
    print(f"Spinless pairing: {encoded.num_qubits} qubits, {len(encoded.terms)} Pauli terms; E/t={energies[0]:.10f}")
    print(f"⟨c₀c₁⟩={pair_expectation:.6f}; Var(N)={variance:.6f}; parity conserved")
    print("Full matrix, state layout and observables checked against independent occupation-bit CAR.")
    finish(args, data, started, plot)


if __name__ == "__main__":
    main()
