# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/tests/test_abinitio.py. Changes: pyalps result
# names (eigenvalues, eigenvectors, problem.states); the Hermiticity check
# uses the fock_matrix fixture; the engine-alias test is gone with the
# selector, replaced by a sector check of several states.
"""Hamiltonian construction and the determinant/parity bridge.

The parity bridge is the load-bearing convention in this package: ALPS
orders fermion modes interleaved (dn before up within each site), while
qdk-chemistry, MACIS and PySCF order them blocked (all alpha, then all
beta). A sign error here is silent and corrupts every CI vector, so it is
pinned both by unit cases and against PySCF FCI.
"""

from __future__ import annotations

import numpy as np
import pytest

from qdk_chemistry_alps.abinitio import (
    alps_state_to_occupations,
    build_abinitio_hamiltonian,
    interleaved_to_blocked_parity,
    solve_ground_state,
)

# spinful_fermion() local basis: 0=|0>, 1=|up>, 2=|dn>, 3=|up,dn>.
LOCAL_STATES = [
    ([0, 0], ([0, 0], [0, 0])),
    ([1, 0], ([1, 0], [0, 0])),
    ([2, 0], ([0, 0], [1, 0])),
    ([3, 0], ([1, 0], [1, 0])),
    ([1, 2], ([1, 0], [0, 1])),
    ([3, 3], ([1, 1], [1, 1])),
]


@pytest.mark.parametrize(("row", "expected"), LOCAL_STATES)
def test_local_states_decode_to_occupations(row, expected):
    alpha, beta = alps_state_to_occupations(np.array(row))
    assert alpha.tolist() == expected[0]
    assert beta.tolist() == expected[1]


def test_local_basis_matches_the_decoding():
    """The decoding tables follow the pyalps site, not a remembered layout."""
    from pyalps import models

    site = models.spinful_fermion()
    n_up = np.diag(np.asarray(site.op("n_up").matrix).real)
    n_dn = np.diag(np.asarray(site.op("n_dn").matrix).real)
    for local in range(4):
        alpha, beta = alps_state_to_occupations(np.array([local]))
        assert (alpha[0], beta[0]) == (n_up[local], n_dn[local])
    # |up,dn> = c+_dn c+_up |0>: the dn mode is created last.
    assert np.asarray(site.op("c_dag_up").matrix)[3, 2] == -1
    assert np.asarray(site.op("c_dag_dn").matrix)[3, 1] == 1


def _brute_force_parity(alpha, beta) -> int:
    """Independent reference: parity of the permutation, counted directly.

    Writes out the occupied-mode string in ALPS order (per site, dn before
    up) and in blocked order (all alpha ascending, then all beta), then
    counts inversions between them.
    """
    interleaved = []
    for i in range(len(alpha)):
        if beta[i]:
            interleaved.append(("d", i))
        if alpha[i]:
            interleaved.append(("u", i))
    blocked = [("u", i) for i in range(len(alpha)) if alpha[i]]
    blocked += [("d", i) for i in range(len(beta)) if beta[i]]

    permutation = [interleaved.index(mode) for mode in blocked]
    inversions = sum(
        1
        for a in range(len(permutation))
        for b in range(a + 1, len(permutation))
        if permutation[a] > permutation[b]
    )
    return -1 if inversions % 2 else 1


@pytest.mark.parametrize(
    ("alpha", "beta", "sign"),
    [
        ([0, 0], [0, 0], +1),   # vacuum: nothing to reorder
        ([1, 0], [0, 0], +1),   # alpha only: already blocked
        ([0, 0], [1, 0], +1),   # beta only: already blocked
        ([1, 0], [1, 0], -1),   # doubly occupied: one dn-before-up swap
        ([0, 1], [1, 0], -1),   # d(0) u(1) -> u(1) d(0): one crossing
        ([1, 0], [0, 1], +1),   # u(0) d(1): already blocked
        ([1, 1], [1, 1], -1),   # two doubly occupied orbitals: three swaps
    ],
)
def test_parity_bridge_signs(alpha, beta, sign):
    assert interleaved_to_blocked_parity(np.array(alpha), np.array(beta)) == sign
    assert _brute_force_parity(alpha, beta) == sign


def test_parity_matches_a_brute_force_permutation_count():
    """Cross-check the closed-form suffix-sum against direct inversion counting."""
    rng = np.random.default_rng(0)
    for _ in range(200):
        n = int(rng.integers(1, 7))
        alpha = rng.integers(0, 2, size=n)
        beta = rng.integers(0, 2, size=n)
        assert interleaved_to_blocked_parity(alpha, beta) == _brute_force_parity(
            alpha.tolist(), beta.tolist()
        )


# --- Hamiltonian construction ---------------------------------------------

def test_builder_validates_shapes():
    with pytest.raises(ValueError, match="h1 must be square"):
        build_abinitio_hamiltonian(np.zeros((2, 3)), np.zeros((2, 2, 2, 2)))
    with pytest.raises(ValueError, match=r"eri must be"):
        build_abinitio_hamiltonian(np.zeros((2, 2)), np.zeros((3, 3, 3, 3)))


def test_builder_drops_terms_below_the_cutoff(integrals):
    """A large cutoff must produce a strictly smaller operator."""
    h1, eri = integrals(3)
    dense = build_abinitio_hamiltonian(h1, eri, tol=0.0)
    sparse = build_abinitio_hamiltonian(h1, eri, tol=1e6)
    assert len(sparse.expand().terms) < len(dense.expand().terms)


def test_hamiltonian_is_hermitian(integrals, fock_matrix):
    h1, eri = integrals(3)
    matrix = fock_matrix(build_abinitio_hamiltonian(h1, eri))
    assert matrix.shape == (64, 64)
    assert np.max(np.abs(matrix - matrix.conj().T)) < 1e-12


# --- against PySCF ---------------------------------------------------------

def _mo_integrals(mol):
    from pyscf import ao2mo, scf

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    c = mf.mo_coeff
    n = c.shape[1]
    h1 = c.T @ mf.get_hcore() @ c
    eri = ao2mo.restore(1, ao2mo.kernel(mol, c), n)
    return h1, eri, float(mol.energy_nuc())


@pytest.mark.pyscf
@pytest.mark.parametrize(
    ("name", "atom", "basis", "n_alpha", "n_beta"),
    [
        ("H2", "H 0 0 0; H 0 0 0.74", "sto-3g", 1, 1),
        ("H4", "H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0", "sto-3g", 2, 2),
    ],
)
def test_ground_state_energy_matches_pyscf_fci(name, atom, basis, n_alpha, n_beta):
    from pyscf import fci, gto

    mol = gto.M(atom=atom, basis=basis, verbose=0)
    h1, eri, ecore = _mo_integrals(mol)

    H = build_abinitio_hamiltonian(h1, eri)
    energy = float(solve_ground_state(H, n_alpha, n_beta).eigenvalues[0]) + ecore
    reference = fci.direct_spin1.FCI().kernel(
        h1, eri, h1.shape[0], (n_alpha, n_beta)
    )[0] + ecore

    assert abs(energy - reference) < 1e-8, name


@pytest.mark.pyscf
def test_ci_vector_is_sign_exact_against_pyscf():
    """The parity bridge, end to end: coefficients agree including sign."""
    from pyscf import fci, gto

    mol = gto.M(
        atom="H 0 0 0; H 0 0 1.0; H 0 0 2.0; H 0 0 3.0", basis="sto-3g", verbose=0
    )
    h1, eri, _ = _mo_integrals(mol)
    n, n_alpha, n_beta = h1.shape[0], 2, 2

    result = solve_ground_state(build_abinitio_hamiltonian(h1, eri), n_alpha, n_beta)
    vec = np.asarray(result.eigenvectors)[:, 0]

    _, ci_ref = fci.direct_spin1.FCI().kernel(h1, eri, n, (n_alpha, n_beta))
    strings_a = fci.cistring.make_strings(range(n), n_alpha)
    address = {int(s): i for i, s in enumerate(strings_a)}

    ci = np.zeros_like(ci_ref)
    for i, row in enumerate(np.asarray(result.problem.states)):
        alpha, beta = alps_state_to_occupations(row)
        key_a = address[int(sum(1 << p for p in np.nonzero(alpha)[0]))]
        key_b = address[int(sum(1 << p for p in np.nonzero(beta)[0]))]
        ci[key_a, key_b] = vec[i] * interleaved_to_blocked_parity(alpha, beta)

    # Fix the arbitrary global phase on the largest reference amplitude.
    pivot = np.unravel_index(np.argmax(np.abs(ci_ref)), ci_ref.shape)
    ci *= np.sign(ci_ref[pivot]) * np.sign(ci[pivot])
    assert np.max(np.abs(ci - ci_ref)) < 1e-6


def test_ground_state_sector_holds_every_requested_state(integrals):
    h1, eri = integrals(3)
    result = solve_ground_state(build_abinitio_hamiltonian(h1, eri), 2, 1, k=2)
    occupations = [alps_state_to_occupations(row) for row in result.problem.states]
    assert all(a.sum() == 2 and b.sum() == 1 for a, b in occupations)
    # 3 alpha placements x 3 beta placements; two ascending energies.
    assert len(occupations) == 9
    assert np.asarray(result.eigenvectors).shape == (9, 2)
    assert result.eigenvalues[0] <= result.eigenvalues[1]
