# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/src/qdk_chemistry_alps/abinitio.py. Changes: the
# Hamiltonian is a pyalps.models.Hamiltonian on spinful_fermion() sites and
# zero integrals are always skipped; solve_ground_state returns the
# pyalps.ed eigenstates (eigenvalues, eigenvectors, problem.states) and has
# no engine selector, whose names all meant the one native solver. The
# determinant helpers are unchanged, except that interleaved_to_blocked_parity
# also takes stacked rows (written here); pyalps's spinful_fermion() has
# alpsx's local basis and matrices. The term loops now run natively in
# pyalps.models.ab_initio (same terms, order and coefficients; written here).
"""Ab initio Hamiltonians on the pyalps operator algebra.

Builds a :class:`pyalps.models.Hamiltonian` from molecular integrals
(chemist's notation) and solves it with :mod:`pyalps.ed`.

Conventions
-----------
- ``h1[p, q]``  : one-electron integrals h_pq.
- ``eri[p, q, r, s]`` : two-electron integrals (pq|rs), chemist's notation,
  8-fold symmetric for real orbitals.
- Second-quantized Hamiltonian::

      H = sum_{pq,s} h_pq c^+_{ps} c_{qs}
        + 1/2 sum_{pqrs,st} (pq|rs) c^+_{ps} c^+_{rt} c_{st} c_{qs}
        + E_core

The determinant helpers at the bottom of this module
(:func:`alps_state_to_occupations`, :func:`interleaved_to_blocked_parity`)
are part of the supported surface: they are the bridge between ALPS's
interleaved mode ordering and the blocked alpha-then-beta ordering that
qdk-chemistry, MACIS and PySCF use, and they are what to call when
checking a sign convention independently.
"""

from __future__ import annotations

import numpy as np
from pyalps import ed, models

__all__ = [
    "alps_state_to_occupations",
    "build_abinitio_hamiltonian",
    "interleaved_to_blocked_parity",
    "solve_ground_state",
]

def build_abinitio_hamiltonian(
    h1: np.ndarray,
    eri: np.ndarray,
    *,
    tol: float = 1e-12,
) -> models.Hamiltonian:
    """Return a pyalps Hamiltonian for the given molecular integrals.

    One site per spatial orbital (``spinful_fermion()``, dimension 4) on an
    open chain, under the strict fermion convention with ``N`` and ``Sz``
    conserved. Every integral becomes an explicit operator term, so no
    lattice geometry is implied or used. Integrals that are zero or smaller
    in magnitude than ``tol`` are dropped.
    """
    h1 = np.asarray(h1, dtype=float)
    n = h1.shape[0] if h1.ndim else 0
    if h1.shape != (n, n):
        raise ValueError(f"h1 must be square, got {h1.shape}")
    eri = np.asarray(eri, dtype=float)
    if eri.shape != (n, n, n, n):
        raise ValueError(f"eri must be ({n},)*4, got {eri.shape}")
    # One-body: sum_{pq,s} h_pq c^+_{ps} c_{qs}; two-body: 1/2 sum_{pqrs,st}
    # (pq|rs) c^+_{ps} c^+_{rt} c_{st} c_{qs}, both built natively.
    return models.ab_initio(h1, eri, tol=tol)


def solve_ground_state(H: models.Hamiltonian, n_alpha: int, n_beta: int, *, k: int = 1):
    """Diagonalize H in the (N, Sz) sector of ``n_alpha`` and ``n_beta``.

    Returns the :mod:`pyalps.ed` eigenstates: ascending ``eigenvalues``,
    ``eigenvectors`` with one column per state, and ``problem``, whose
    ``states`` list the per-site local states of each sector basis row.
    """
    sector = {"N": float(n_alpha + n_beta), "Sz": 0.5 * (n_alpha - n_beta)}
    return ed.prepare_product(H, sector=sector).eigenstates(count=k)


# ---------------------------------------------------------------------------
# Determinant bookkeeping: ALPS basis rows <-> (alpha, beta) occupations
# ---------------------------------------------------------------------------
# Local basis of spinful_fermion(): 0=|0>, 1=|up>, 2=|dn>, 3=|up,dn>.
# The matrices (c_dag_up[3,2] = -1, c_dag_dn[3,1] = +1) realize
# |up,dn> = c^+_dn c^+_up |0> — the dn mode comes FIRST within each site.
# The global Fock ordering is therefore mode-interleaved as
#   d(0), u(0), d(1), u(1), ..., d(n-1), u(n-1)
# (sites ascending: the strict convention places Jordan-Wigner strings on
# the sites below each odd-parity factor).
# PySCF orders all alpha modes before all beta modes, ascending orbital.

_LOCAL_UP = np.array([0, 1, 0, 1], dtype=np.int8)   # n_up per local state
_LOCAL_DN = np.array([0, 0, 1, 1], dtype=np.int8)   # n_dn per local state


def alps_state_to_occupations(
    state_row: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Map one ALPS basis row (per-site local states) to (alpha_occ, beta_occ)."""
    local = np.asarray(state_row, dtype=np.int64)
    return _LOCAL_UP[local].astype(np.int8), _LOCAL_DN[local].astype(np.int8)


def interleaved_to_blocked_parity(alpha_occ: np.ndarray, beta_occ: np.ndarray) -> int | np.ndarray:
    """Sign relating the interleaved mode ordering to PySCF's blocked ordering.

    One determinant (1-D occupations) gives an int; stacked determinants
    (occupations of shape ``(..., n)``) give an int64 array of signs.

    Reordering the creation-operator string from
      d(0) u(0) d(1) u(1) ...   (ALPS/JW order, occupied modes only)
    to
      u(0) u(1) ... d(0) d(1) ...  (PySCF alpha-block-then-beta-block)
    costs one transposition per (occupied beta on orbital p, occupied alpha
    on orbital q >= p) pair: cross-site pairs with q > p, plus one pair per
    doubly-occupied orbital from the dn-before-up site-local convention.
    """
    alpha_occ = np.asarray(alpha_occ, dtype=np.int64)
    beta_occ = np.asarray(beta_occ, dtype=np.int64)
    # inclusive suffix count of alpha electrons at orbitals q >= p
    alpha_from = np.cumsum(alpha_occ[..., ::-1], axis=-1)[..., ::-1]
    swaps = (beta_occ * alpha_from).sum(axis=-1)
    if np.ndim(swaps) == 0:
        return -1 if int(swaps) % 2 else 1
    return 1 - 2 * (swaps % 2)
