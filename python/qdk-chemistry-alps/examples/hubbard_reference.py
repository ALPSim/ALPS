# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/hubbard_reference.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 2fb41daf147809f23ed802194868acead995044e240a7ee7471cc4abfc05318f.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""QDK → ALPS: Hubbard energies, occupations and spin correlations."""

import time
from itertools import combinations

import numpy as np
from demo_support import finish, fock_matrix, parser
from qdk_chemistry.algorithms import create
from qdk_chemistry.data import LatticeGraph
from qdk_chemistry.utils.model_hamiltonians import create_hubbard_hamiltonian

import qdk_chemistry_alps  # noqa: F401 — register ALPS calculators


def reference(n, interaction):
    """Independent fixed-sector Fock reference; 400 determinants at n=6."""
    strings = [sum(1 << i for i in occupied) for occupied in combinations(range(n), n // 2)]
    states = [up | (down << n) for up in strings for down in strings]
    hopping = [(-1., [(offset + i, True), (offset + j, False)])
               for offset in (0, n) for i in range(n) for j in range(n) if abs(i - j) == 1]
    h = fock_matrix(states, hopping)
    occupation = np.array([[(state >> i) & 1 for i in range(2 * n)] for state in states])
    h += np.diag(interaction * np.sum(occupation[:, :n] * occupation[:, n:], axis=1))
    energies, vectors = np.linalg.eigh(h)
    psi = vectors[:, 0]
    sz = .5 * (occupation[:, :n] - occupation[:, n:])
    correlations = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            transverse = fock_matrix(states, [
                (.5, [(i, True), (n + i, False), (n + j, True), (j, False)]),
                (.5, [(n + i, True), (i, False), (j, True), (n + j, False)]),
            ])
            correlations[i, j] = np.vdot(psi, transverse @ psi).real + np.dot(abs(psi)**2, sz[:, i] * sz[:, j])
    density = abs(psi)**2 @ (occupation[:, :n] + occupation[:, n:])
    return energies[0], density, correlations


def plot(plt, data):
    figure, axes = plt.subplots(1, 2, figsize=(9, 3.5), layout="constrained")
    axes[0].plot(range(6), data["occupations"], "o-", color="#087e8b", label="ALPS DMRG")
    axes[0].plot(range(6), data["reference_occupations"], "x", color="#c44e52", label="Independent exact")
    axes[0].set(xlabel="Site", ylabel="⟨nᵢ⟩", ylim=(0, 2), title="Half-filled Hubbard chain")
    axes[0].legend()
    image = axes[1].imshow(data["spin_correlations"], cmap="RdBu_r", vmin=-.75, vmax=.75)
    axes[1].set(xlabel="Site j", ylabel="Site i", title="DMRG spin correlations ⟨Sᵢ·Sⱼ⟩")
    figure.colorbar(image, ax=axes[1])
    return figure


def main():
    args = parser(__doc__).parse_args()
    started = time.perf_counter()
    h = create_hubbard_hamiltonian(LatticeGraph.chain(6), epsilon=0., t=1., U=4.)
    ed = create("multi_configuration_calculator", "alps_ed")
    energy_ed, _ = ed.run(h, 3, 3)
    dmrg = create("multi_configuration_calculator", "alps_classic_dmrg")
    for key, value in {"max_bond_dim": 128, "sweeps": 12, 
                       "truncation_tol": 0.0}.items():
        dmrg.settings().set(key, value)
    energy_dmrg, _ = dmrg.run(h, 3, 3)
    occupations = np.asarray(dmrg.measure_local("n")).real
    spin = np.asarray(dmrg.spin_spin_correlations()).real
    density = np.asarray(dmrg.measure_correlations("n")).real
    dm1, dm2 = map(np.asarray, dmrg.make_rdms())
    index = np.arange(6)
    density_from_rdm = dm2[index, index][:, index, index] + np.diag(np.diag(dm1))
    exact, exact_density, exact_spin = reference(6, 4.)
    np.testing.assert_allclose([energy_ed, energy_dmrg], exact, atol=1e-8, rtol=0)
    np.testing.assert_allclose(occupations, exact_density, atol=1e-8, rtol=0)
    np.testing.assert_allclose(spin, exact_spin, atol=1e-8, rtol=0)
    np.testing.assert_allclose(occupations, np.diag(dm1), atol=1e-9, rtol=0)
    np.testing.assert_allclose(density, density_from_rdm, atol=1e-9, rtol=0)
    assert abs(occupations.sum() - 6) < 1e-9
    print(f"Hubbard, 6 sites, open, N↑=N↓=3, U/t=4\nED   E/t = {energy_ed:.10f}\nDMRG E/t = {energy_dmrg:.10f}")
    print(f"Independent energy/correlation checks passed; max spin error = {np.max(abs(spin - exact_spin)):.2e}")
    print("Nearest-neighbour ⟨Sᵢ·Sᵢ₊₁⟩:", np.round(np.diag(spin, 1), 5))
    finish(args, {"model": {"sites": 6, "t": 1., "U": 4., "n_alpha": 3, "n_beta": 3, "boundary": "open"},
                  "controls": {"max_bond_dim": 128, "sweeps": 12,  "truncation": 0.0},
                  "energy_ed": energy_ed, "energy_dmrg": energy_dmrg, "energy_reference": exact,
                  "occupations": occupations, "reference_occupations": exact_density,
                  "spin_correlations": spin, "reference_spin_correlations": exact_spin}, started, plot)


if __name__ == "__main__":
    main()
