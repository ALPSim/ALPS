# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/src/qdk_chemistry_alps/interop.py. Changes:
# pyalps names throughout (models.lattice.UnitCell and Lattice.build for the
# graph lattice, models.spin() and models.located for spin operators,
# Factor.site/op/basis and models.Statistics for expanded terms,
# models.PauliOperator from models.export_pauli); hamiltonian_to_integrals
# decomposes each local factor matrix into the site's own ladder operators
# instead of looking its decomposition up by operator name, so derived names
# (c_up_dag, ...) and other spinful sites with c_dag_up and c_dag_dn need no
# table. The normal-ordering engine, the integral extraction and its
# validation, the qubit mappings and the state placement are unchanged.
"""Thin adapters between qdk-chemistry's data containers and the pyalps
lattice / operator-algebra representations.

qdk-chemistry's lattice surface is a pure adjacency graph
(:class:`~qdk_chemistry.data.LatticeGraph`); its spin models travel as
Pauli-string :class:`~qdk_chemistry.data.QubitOperator` (little-endian:
the character at position ``n-1-i`` acts on qubit ``i``); its fermionic
models travel as integral Hamiltonians (h1/eri, chemist notation). The
pyalps side has geometric :class:`pyalps.models.lattice.Lattice` objects
and a symbolic operator algebra. These adapters translate at the boundary
so that neither package's model factories need duplicating:

- :func:`lattice_graph_to_alps` / :func:`alps_lattice_to_lattice_graph`
- :func:`qubit_operator_to_hamiltonian` / :func:`hamiltonian_to_qubit_operator`
  (spin-1/2 / dim-2 bosonic sites)
- :func:`pauli_operator_to_qubit_operator` — wrap a native
  :func:`pyalps.models.export_pauli` expansion (finite spins and bosons,
  strict fermions) as a QubitOperator
- :func:`hamiltonian_to_integrals` — extract (h1, eri, ecore) from a
  particle-conserving strict-fermion pyalps Hamiltonian, with a
  tolerance-controlled round-trip validation (so qdk's own solvers, e.g.
  MACIS, can consume pyalps-built models)
- :func:`fermionic_hamiltonian_to_qubit_operator` — the qubit-operator export
  route for restricted fermionic models: pyalps Hamiltonian -> integrals ->
  qdk's Majorana-loop mapping engine -> QubitOperator, for supported
  untapered encodings (Jordan-Wigner, Bravyi-Kitaev, parity, ...)
- :func:`ci_vector_to_qubit_statevector` — place an ALPS CI vector into the
  blocked Jordan–Wigner register of a chemistry-route export, for
  expectation-value comparisons on an ordinary computer; this does not
  synthesize a state-preparation circuit
- :func:`rcm_ordering` — bandwidth-reducing DMRG site ordering from a
  LatticeGraph
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # imported lazily at runtime inside each function
    from pyalps.models import Hamiltonian, PauliOperator
    from pyalps.models.lattice import Lattice
    from qdk_chemistry.data import LatticeGraph, MajoranaMapping, QubitOperator

__all__ = [
    "lattice_graph_to_alps",
    "alps_lattice_to_lattice_graph",
    "qubit_operator_to_hamiltonian",
    "hamiltonian_to_qubit_operator",
    "pauli_operator_to_qubit_operator",
    "hamiltonian_to_integrals",
    "fermionic_hamiltonian_to_qubit_operator",
    "ci_vector_to_qubit_statevector",
    "rcm_ordering",
]


# ---------------------------------------------------------------------------
# Lattices
# ---------------------------------------------------------------------------

def lattice_graph_to_alps(graph: LatticeGraph) -> tuple[Lattice, list[float]]:
    """qdk LatticeGraph -> pyalps custom Lattice.

    Every distinct edge weight becomes its own bond type (ascending weight
    order), so weighted graphs map onto pyalps factories that accept
    per-bond-type couplings. Returns ``(lattice, type_weights)`` where
    ``type_weights[t]`` is the edge weight carried by bond type ``t``.
    A LatticeGraph has no geometry, so the sites form one open cell with
    1-D coordinates equal to their index.
    """
    from pyalps.models import lattice as lat

    A = np.asarray(graph.adjacency_matrix(), dtype=float)
    n = A.shape[0] if A.ndim == 2 else 0
    if n == 0 or A.shape != (n, n) or not np.allclose(A, A.T):
        raise ValueError("LatticeGraph adjacency must be square and symmetric")

    src, tgt = np.nonzero(np.triu(A, k=1))
    weights = A[src, tgt]
    type_weights = sorted(set(float(w) for w in weights))
    type_of = {w: t for t, w in enumerate(type_weights)}

    cell = lat.UnitCell(1)
    for i in range(n):
        # Fractional coordinates of a cell of length n: site i sits at x = i.
        cell.add_vertex([i / n])
    for s, t, w in zip(src, tgt, weights, strict=True):
        cell.add_edge(int(s), int(t), type_of[float(w)])
    lattice = lat.Lattice.build(cell, [[float(n)]], [1], "open")
    return lattice, type_weights


def alps_lattice_to_lattice_graph(
    lattice: Lattice,
    type_weights: list[float] | None = None,
) -> LatticeGraph:
    """pyalps Lattice -> qdk LatticeGraph (adjacency only; geometry and
    periodic wrapping information are not expressible in a LatticeGraph).

    ``type_weights`` optionally maps bond types to edge weights
    (default: every edge weight 1.0).
    """
    from qdk_chemistry.data import LatticeGraph

    n = lattice.num_sites
    A = np.zeros((n, n))
    for b in lattice.bonds:
        w = 1.0 if type_weights is None else float(type_weights[b.bond_type])
        A[b.source, b.target] = w
        A[b.target, b.source] = w
    return LatticeGraph.from_dense_matrix(A)


def rcm_ordering(graph: LatticeGraph) -> np.ndarray:
    """Reverse Cuthill-McKee site ordering for a LatticeGraph.

    Returns a permutation ``order`` such that relabelling site
    ``order[k] -> k`` minimizes the adjacency bandwidth — the standard
    heuristic for a good 1-D DMRG site ordering of a 2-D/irregular graph.
    Apply with ``adj[order][:, order]`` (and permute integrals the same
    way for fermionic models).
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import reverse_cuthill_mckee

    A = csr_matrix(np.asarray(graph.adjacency_matrix(), dtype=float))
    return np.asarray(reverse_cuthill_mckee(A, symmetric_mode=True), dtype=int)


# ---------------------------------------------------------------------------
# Spin models: QubitOperator <-> pyalps Hamiltonian
# ---------------------------------------------------------------------------

def qubit_operator_to_hamiltonian(
    qubit_op: QubitOperator,
    lattice: Lattice | None = None,
    conserved_qns: tuple[str, ...] = (),
) -> Hamiltonian:
    """qdk QubitOperator -> pyalps Hamiltonian on spin-1/2 sites.

    qdk Pauli strings are little-endian (rightmost character = qubit 0);
    qubit ``i`` maps to site ``i``. Pauli matrices are ``2*S`` on the
    ``models.spin()`` site, whose local state 0 is spin up (Pauli ``|0>``).
    An all-identity string becomes an explicit ``c * I(0)`` term so constant
    offsets survive the round trip.

    ``lattice`` defaults to an open chain of ``num_qubits`` sites (the
    geometry is irrelevant for global operator terms); pass
    ``conserved_qns=("Sz",)`` only if every string conserves Sz.

    State-index caveat: qdk dense matrices put qubit 0 in the LEAST
    significant bit, pyalps product bases put site 0 in the MOST
    significant position — the operators are identical, but explicit
    matrices/state vectors relate by bit reversal.
    """
    from pyalps import models

    d = qubit_op.to_dict()
    strings = list(d["pauli_strings"])
    cre = np.asarray(d["coefficients"]["real"], dtype=float)
    cim = np.asarray(d["coefficients"]["imag"], dtype=float)

    n = int(qubit_op.num_qubits)
    lat = lattice if lattice is not None else models.lattice.chain(n, bc="open")
    if lat.num_sites != n:
        raise ValueError(
            f"lattice has {lat.num_sites} sites but the operator acts on "
            f"{n} qubits"
        )
    s = models.spin()
    H = models.Hamiltonian(lat, s, conserved_quantum_numbers=list(conserved_qns))

    def pauli(name: str):
        return lambda i: 2.0 * models.located(s, name, i)

    factor = {"X": pauli("Sx"), "Y": pauli("Sy"), "Z": pauli("Sz")}

    for string, re_c, im_c in zip(strings, cre, cim, strict=True):
        if len(string) != n:
            raise ValueError(f"Pauli string {string!r} has wrong length")
        coeff = complex(re_c, im_c)
        if coeff == 0:
            continue
        if abs(coeff.imag) < 1e-14 * max(1.0, abs(coeff.real)):
            coeff = coeff.real
        term = None
        for pos, ch in enumerate(string):
            if ch == "I":
                continue
            site = n - 1 - pos  # little-endian string -> site index
            f = factor.get(ch)
            if f is None:
                raise ValueError(f"unknown Pauli character {ch!r} in {string!r}")
            piece = f(site)
            term = piece if term is None else term * piece
        if term is None:  # all-identity: constant offset
            term = models.located(s, "I", 0)
        H.add(coeff * term)
    return H


_PAULI_MATS = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_operator_to_qubit_operator(operator: PauliOperator) -> QubitOperator:
    """Convert a native pyalps Pauli expansion into a QDK container.

    Build ``operator`` with :func:`pyalps.models.export_pauli`. This direct
    route supports binary/Gray finite spin and boson spaces, spinless pairing
    and complex/spin-mixing spinful models without chemistry integrals. Keep
    the pyalps result alongside the QDK operator: its ``encoding``,
    ``fermion_mode_qubits`` and ``state_index`` define the
    state/observable layout. QDK's blocked/interleaved chemistry metadata
    cannot describe this local-basis convention, so this adapter
    deliberately leaves those fields unset. Do not use
    ``ci_vector_to_qubit_statevector`` with this route.

    An explicit zero identity preserves the register size for a zero
    operator. Non-Hermitian observable coefficients are retained.
    """
    from pyalps.models import PauliOperator
    from qdk_chemistry.data import QubitOperator

    if not isinstance(operator, PauliOperator):
        raise TypeError("expected pyalps.models.PauliOperator")
    terms = operator.terms
    if not terms:
        return QubitOperator(["I" * operator.num_qubits], np.zeros(1))
    return QubitOperator(list(terms), np.asarray(list(terms.values()), dtype=complex))


def hamiltonian_to_qubit_operator(
    H: Hamiltonian, tol: float = 1e-12
) -> QubitOperator:
    """pyalps Hamiltonian on dim-2 bosonic sites -> qdk QubitOperator.

    Each term's per-site operator product is collapsed to a 2x2 matrix and
    decomposed in the Pauli basis (c_P = tr(P M) / 2); the cartesian
    expansion over the term's support accumulates into Pauli strings
    (little-endian, matching qdk's convention). Fermionic sites are
    rejected — use :func:`fermionic_hamiltonian_to_qubit_operator`,
    which chooses a fermion-to-qubit encoding.
    """
    from pyalps import models
    from qdk_chemistry.data import QubitOperator

    n = H.num_sites
    acc: dict[str, complex] = {}

    for term in H.expand().terms:
        coeff = complex(term.effective_coeff)
        site_mats: dict[int, np.ndarray] = {}
        for factor in term.factors:
            basis = factor.basis
            if basis.dim != 2 or basis.statistics != models.Statistics.bosonic:
                raise ValueError(
                    "hamiltonian_to_qubit_operator supports dim-2 bosonic "
                    f"(spin-1/2) sites only; site {factor.site} has dimension "
                    f"{basis.dim} and {basis.statistics} statistics"
                )
            m = np.asarray(factor.op.matrix, dtype=complex)
            site_mats[factor.site] = site_mats.get(factor.site, np.eye(2)) @ m

        # Pauli decomposition per involved site.
        decomps = []
        for site_idx, m in sorted(site_mats.items()):
            cs = [(p, np.trace(_PAULI_MATS[p].conj().T @ m) / 2.0)
                  for p in ("I", "X", "Y", "Z")]
            decomps.append((site_idx, [(p, c) for p, c in cs if abs(c) > tol]))

        assignments = [{}]
        weights = [coeff]
        for site_idx, cs in decomps:
            new_assign, new_w = [], []
            for a, w in zip(assignments, weights, strict=True):
                for p, c in cs:
                    na = dict(a)
                    if p != "I":
                        na[site_idx] = p
                    new_assign.append(na)
                    new_w.append(w * c)
            assignments, weights = new_assign, new_w

        for a, w in zip(assignments, weights, strict=True):
            if abs(w) <= tol:
                continue
            chars = ["I"] * n
            for site_idx, p in a.items():
                chars[n - 1 - site_idx] = p  # little-endian
            key = "".join(chars)
            acc[key] = acc.get(key, 0.0) + w

    strings = [k for k, v in acc.items() if abs(v) > tol]
    coeffs = np.array([acc[k] for k in strings])
    if np.abs(coeffs.imag).max(initial=0.0) < tol:
        coeffs = coeffs.real
    return QubitOperator(strings, coeffs)


# ---------------------------------------------------------------------------
# Fermionic models: pyalps Hamiltonian -> (h1, eri, ecore)
# ---------------------------------------------------------------------------
# A spinful site's operators span the 16 products O_up O_dn with
# O in (1, c+, c, n = c+ c) of each spin mode. Daggers of each choice:
_ONE_MODE = ((), (True,), (False,), (True, False))
_LADDER_CACHE: dict = {}


def _ladder_decomposer(site_idx, basis):
    """Ladder decomposition of a spinful site's local matrices.

    Returns ``(solve, sequences)``: ``solve(matrix)`` gives the 16
    coefficients of ``matrix`` on the products ``O_up O_dn``, built from the
    site's own ``c_dag_up`` and ``c_dag_dn`` (so the site's local state order
    and intra-site sign convention are respected), and ``sequences[k]`` the
    (spin, dagger) ladder string of product ``k`` (spin 0 = up, 1 = dn).
    """
    if (
        basis.dim != 4
        or basis.statistics != _fermionic()
        or not (basis.has_op("c_dag_up") and basis.has_op("c_dag_dn"))
    ):
        raise ValueError(
            "hamiltonian_to_integrals expects spinful-fermion sites with "
            "c_dag_up and c_dag_dn (pyalps.models.spinful_fermion()); "
            f"site {site_idx} has dimension {basis.dim} and "
            f"{basis.statistics} statistics"
        )
    up = np.asarray(basis.op("c_dag_up").matrix, dtype=complex)
    dn = np.asarray(basis.op("c_dag_dn").matrix, dtype=complex)
    key = (up.tobytes(), dn.tobytes())
    cached = _LADDER_CACHE.get(key)
    if cached is not None:
        return cached

    def ops(cdag):
        return (np.eye(4), cdag, cdag.conj().T, cdag @ cdag.conj().T)

    columns, sequences = [], []
    for a, o_up in enumerate(ops(up)):
        for b, o_dn in enumerate(ops(dn)):
            columns.append((o_up @ o_dn).reshape(-1))
            sequences.append([(0, d) for d in _ONE_MODE[a]] + [(1, d) for d in _ONE_MODE[b]])
    products = np.column_stack(columns)
    if np.linalg.cond(products) > 1e8:
        raise ValueError(
            f"c_dag_up and c_dag_dn of site {site_idx} do not generate its "
            "local operators"
        )
    inverse = np.linalg.inv(products)

    def solve(matrix):
        return inverse @ np.asarray(matrix, dtype=complex).reshape(-1)

    _LADDER_CACHE[key] = (solve, sequences)
    return solve, sequences


def _fermionic():
    from pyalps import models

    return models.Statistics.fermionic


def _normal_order(seq, coeff, out):
    """Accumulate `coeff * prod(seq)` into `out` in canonical normal order.

    seq: tuple of (mode, dagger). Canonical form: all daggers first with
    strictly ascending modes, then annihilators with strictly descending
    modes. Uses c_i c+_j = delta_ij - c+_j c_i and anticommutation;
    repeated identical operators annihilate the term.
    """
    stack = [(list(seq), coeff)]
    while stack:
        ops, c = stack.pop()
        for i in range(len(ops) - 1):
            (m1, d1), (m2, d2) = ops[i], ops[i + 1]
            if not d1 and d2:  # annihilator before creator: contract + swap
                if m1 == m2:
                    stack.append((ops[:i] + ops[i + 2:], c))
                stack.append((ops[:i] + [ops[i + 1], ops[i]] + ops[i + 2:], -c))
                break
            if d1 == d2 and m1 == m2:  # c+c+ or cc on the same mode
                break  # term is zero
            swap = (d1 and d2 and m1 > m2) or (not d1 and not d2 and m1 < m2)
            if swap:
                stack.append((ops[:i] + [ops[i + 1], ops[i]] + ops[i + 2:], -c))
                break
        else:
            key = tuple((m, d) for m, d in ops)
            out[key] = out.get(key, 0.0) + c
    return out


def _accumulate_hamiltonian_terms(terms, n_orb):
    """Normal-order a sequence of (coeff, [(mode, dag), ...]) products."""
    acc: dict = {}
    for coeff, seq in terms:
        _normal_order(seq, coeff, acc)
    return {k: v for k, v in acc.items() if abs(v) > 1e-14}


def _chemist_products(h1, eri, tol=1e-14):
    """(h1, eri) -> the flat (coeff, ladder-sequence) products of
    H = sum h_pq c+_ps c_qs + 1/2 sum (pq|rs) c+_ps c+_rt c_st c_qs,
    modes indexed m = 2*orbital + spin (0=up, 1=dn)."""
    n = h1.shape[0]
    out = []
    for p in range(n):
        for q in range(n):
            if abs(h1[p, q]) > tol:
                for s in (0, 1):
                    out.append((h1[p, q],
                                [(2 * p + s, True), (2 * q + s, False)]))
    for p in range(n):
        for q in range(n):
            for r in range(n):
                for s_orb in range(n):
                    v = eri[p, q, r, s_orb]
                    if abs(v) <= tol:
                        continue
                    for s1 in (0, 1):
                        for s2 in (0, 1):
                            out.append((0.5 * v, [
                                (2 * p + s1, True), (2 * r + s2, True),
                                (2 * s_orb + s2, False), (2 * q + s1, False),
                            ]))
    return out


def hamiltonian_to_integrals(
    H: Hamiltonian, tol: float = 1e-10
) -> tuple[np.ndarray, np.ndarray, float]:
    """Extract ``(h1, eri, ecore)`` from a particle-conserving
    strict-fermion pyalps Hamiltonian on spinful-fermion sites.

    Each local factor is decomposed into the site's ladder operators and
    the Hamiltonian's terms are normal-ordered by a small fermionic Wick
    engine; the one-body block must be spin-diagonal and spin-symmetric,
    and the two-body block must be expressible as
    ``1/2 sum (pq|rs) c+_ps c+_rt c_st c_qs`` with a real 8-fold-symmetric
    ``eri`` (the chemist form qdk-chemistry's containers carry). The
    extraction is validated within ``tol``: the chemist form of the extracted
    integrals is re-normal-ordered and compared term by term against the
    input — any Hamiltonian outside this class (pairing terms, spin flips,
    genuine 3-body, spin-asymmetric interactions) raises ``ValueError``.
    """
    n = H.num_sites

    products = []
    for term in H.expand().terms:
        coeff = complex(term.effective_coeff)
        if abs(coeff.imag) > 1e-14 * max(1.0, abs(coeff.real)):
            raise ValueError("complex coefficients are not chemist-representable")
        pieces = [(coeff.real, [])]
        for factor in term.factors:
            solve, sequences = _ladder_decomposer(factor.site, factor.basis)
            weights = solve(factor.op.matrix)
            if np.abs(weights.imag).max() > 1e-14 * max(1.0, np.abs(weights).max()):
                raise ValueError(
                    f"operator {factor.op.name!r} on site {factor.site} has "
                    "complex ladder coefficients, which are not "
                    "chemist-representable"
                )
            elem = [
                (w.real, [(2 * factor.site + spin, dag) for spin, dag in seq])
                for w, seq in zip(weights, sequences, strict=True)
                if abs(w) > 1e-15
            ]
            pieces = [
                (c0 * c1, seq0 + seq1)
                for c0, seq0 in pieces
                for c1, seq1 in elem
            ]
        products.extend(pieces)

    acc = _accumulate_hamiltonian_terms(products, n)

    ecore = float(acc.get((), 0.0))
    h1 = np.zeros((n, n))
    eri = np.zeros((n, n, n, n))
    for key, c in acc.items():
        ndag = sum(1 for _, d in key if d)
        nann = len(key) - ndag
        if ndag != nann:
            raise ValueError(
                f"term {key} does not conserve particle number"
            )
        if len(key) == 0:
            continue
        if len(key) == 2:
            (m1, _), (m2, _) = key
            p, s1 = divmod(m1, 2)
            q, s2 = divmod(m2, 2)
            if s1 != s2:
                raise ValueError(f"one-body term {key} flips spin")
            h1[p, q] += c / 1.0  # accumulated per spin; validated below
        elif len(key) == 4:
            pass  # handled below via the cross-spin block
        else:
            raise ValueError(
                f"term {key} is beyond two-body (rank {len(key)})"
            )

    # h1 was accumulated once per spin channel; the two channels must agree.
    h1 /= 2.0

    # eri from the cross-spin normal-ordered block: the canonical
    # coefficient of c+_(p,up) c+_(r,dn) c_(s,dn) c_(q,up) equals
    # sign * (pq|rs) (both spin assignments of the chemist double sum land
    # on it, cancelling the 1/2). _normal_order on the already
    # normal-ordered probe sequence yields exactly the sorting sign.
    for p in range(n):
        for q in range(n):
            for r in range(n):
                for s_orb in range(n):
                    probe: dict = {}
                    _normal_order(
                        [(2 * p, True), (2 * r + 1, True),
                         (2 * s_orb + 1, False), (2 * q, False)],
                        1.0, probe,
                    )
                    (key, sign), = probe.items()
                    eri[p, q, r, s_orb] = sign * acc.get(key, 0.0)

    # Exact round-trip validation: re-normal-order the chemist form of the
    # extracted integrals and require term-by-term agreement.
    acc_pred = _accumulate_hamiltonian_terms(
        _chemist_products(h1, eri) + ([(ecore, [])] if ecore else []), n
    )
    keys = set(acc) | set(acc_pred)
    scale = max((abs(v) for v in acc.values()), default=1.0)
    worst = max(
        (abs(acc.get(k, 0.0) - acc_pred.get(k, 0.0)) for k in keys),
        default=0.0,
    )
    if worst > tol * scale:
        raise ValueError(
            "Hamiltonian is not representable as spin-restricted chemist "
            f"integrals (worst normal-ordered term mismatch {worst:.3e}); "
            "it may be spin-asymmetric or contain non-chemist two-body "
            "structure"
        )
    return h1, eri, ecore


# ---------------------------------------------------------------------------
# Fermionic models: pyalps Hamiltonian -> integrals -> QubitOperator
# ---------------------------------------------------------------------------

_MAJORANA_ENCODINGS = {
    "jordan-wigner": "jordan_wigner",
    "bravyi-kitaev": "bravyi_kitaev",
    "bravyi-kitaev-tree": "bravyi_kitaev_tree",
    "parity": "parity",
}


def fermionic_hamiltonian_to_qubit_operator(
    H: Hamiltonian,
    encoding: str | MajoranaMapping = "jordan-wigner",
    *,
    threshold: float = 1e-12,
    integral_threshold: float = 1e-12,
) -> QubitOperator:
    """Particle-conserving strict-fermion pyalps Hamiltonian -> qdk
    QubitOperator, via a fermion-to-qubit encoding.

    The restricted-integral export route for fermionic lattice models:
    :func:`hamiltonian_to_integrals` recovers the spin-restricted chemist
    integrals (checked by rebuilding the operator within tolerance), and
    qdk-chemistry's Majorana-loop engine (``majorana_map_hamiltonian`` — the
    same C++ path qdk's own ``QdkQubitMapper`` uses) maps them to Pauli
    terms. The constant part of the Hamiltonian is carried through as an
    identity term, so the QubitOperator's spectrum on the full 2^(2n) qubit
    space equals the fermionic Fock-space spectrum up to numerical error and
    coefficient pruning. Export does not construct or execute a circuit.

    ``encoding`` is either one of the names in ``_MAJORANA_ENCODINGS``
    ("jordan-wigner", "bravyi-kitaev", "bravyi-kitaev-tree", "parity") or
    a prebuilt :class:`~qdk_chemistry.data.MajoranaMapping` with
    ``num_modes == 2 * H.num_sites`` and no tapering metadata (for
    symmetry-conserving-BK with tapering, map the qdk Hamiltonian with
    ``QdkQubitMapper``, which applies the tapering post-pass).

    Conventions (matching ``QdkQubitMapper``): spin-orbitals are BLOCKED
    (modes 0..n-1 = up/alpha, n..2n-1 = down/beta; site ``i`` is orbital
    ``i``); Pauli labels are little-endian (rightmost character = qubit 0).
    ``threshold`` prunes small Pauli coefficients, ``integral_threshold``
    skips small integrals.
    """
    from qdk_chemistry._core.data import (
        majorana_map_hamiltonian,
        sparse_pauli_word_to_label,
    )
    from qdk_chemistry.data import (
        FermionModeOrder,
        MajoranaMapping,
        QubitOperator,
    )

    h1, eri, ecore = hamiltonian_to_integrals(H)
    n = h1.shape[0]

    if isinstance(encoding, MajoranaMapping):
        mapping = encoding
        if mapping.tapering:
            raise ValueError(
                "the mapping carries tapering metadata, which the raw "
                "engine cannot apply; pass mapping.without_tapering() or "
                "use qdk's QdkQubitMapper for tapered encodings"
            )
        if mapping.num_modes != 2 * n:
            raise ValueError(
                f"MajoranaMapping has {mapping.num_modes} modes but the "
                f"Hamiltonian has {2 * n} spin-orbitals (2 x {n} sites)"
            )
    else:
        ctor = _MAJORANA_ENCODINGS.get(encoding)
        if ctor is None:
            raise ValueError(
                f"unknown encoding {encoding!r}; expected one of "
                f"{sorted(_MAJORANA_ENCODINGS)} or a MajoranaMapping"
            )
        mapping = getattr(MajoranaMapping, ctor)(num_modes=2 * n)

    # The pybind engine does no buffer-shape validation (a wrong-size
    # array is a silent out-of-bounds read), so pin shape/dtype/layout
    # here. spin_symmetric=True: hamiltonian_to_integrals only returns
    # spin-restricted integrals, and the engine still requires all five
    # arrays, so the same objects serve every spin channel.
    h1 = np.ascontiguousarray(h1, dtype=np.float64)
    eri = np.ascontiguousarray(eri, dtype=np.float64)
    if h1.shape != (n, n) or eri.shape != (n, n, n, n):
        raise ValueError("integral arrays have inconsistent shapes")

    words, coeffs = majorana_map_hamiltonian(
        mapping, float(ecore), h1, h1, eri, eri, eri,
        n, True, float(threshold), float(integral_threshold),
    )
    labels = [sparse_pauli_word_to_label(w, mapping.num_qubits) for w in words]
    coeffs = np.asarray(coeffs)
    if coeffs.size and np.abs(coeffs.imag).max() <= threshold:
        coeffs = coeffs.real
    if not labels:  # H below threshold everywhere: explicit zero operator
        labels, coeffs = ["I" * mapping.num_qubits], np.zeros(1)
    return QubitOperator(
        labels, coeffs,
        encoding=mapping.base_encoding,
        fermion_mode_order=FermionModeOrder.BLOCKED,
    )


def ci_vector_to_qubit_statevector(
    basis_states,
    coefficients: np.ndarray,
    n_orbitals: int,
) -> np.ndarray:
    """ALPS CI vector -> amplitudes over the qubit register.

    For the default Jordan–Wigner chemistry export, this places a state
    in ALPS's determinant basis into the corresponding blocked register
    of ``2 * n_orbitals`` qubits. It does not apply the occupation/parity
    transforms needed for Bravyi–Kitaev or parity encoding, and it does not
    implement the separate native-export layout. Keep the state mapping
    matched to the operator encoding when comparing expectation values.

    Two conversions happen here, and both are silent when wrong:

    * **Mode ordering.** ALPS orders fermion modes interleaved; the qubit
      operator assumes blocked. Coefficients pick up the sign returned by
      :func:`~qdk_chemistry_alps.abinitio.interleaved_to_blocked_parity`.
    * **Qubit indexing.** Amplitudes are indexed little-endian — qubit ``q``
      is bit ``q`` of the array index — matching qdk-chemistry, whereas ALPS
      puts site 0 in the most significant position.

    ``basis_states`` are the per-site local states of the sector basis (the
    ``problem.states`` of a :func:`~qdk_chemistry_alps.abinitio.solve_ground_state`
    result) and ``coefficients`` their amplitudes. Returns a complex vector
    of length ``2 ** (2 * n_orbitals)``, zero outside the sector.
    """
    from .abinitio import alps_state_to_occupations, interleaved_to_blocked_parity

    coefficients = np.asarray(coefficients)
    basis_states = np.asarray(basis_states)
    if basis_states.shape[0] != coefficients.shape[0]:
        raise ValueError(
            f"got {basis_states.shape[0]} basis states but "
            f"{coefficients.shape[0]} coefficients"
        )
    if basis_states.shape[1] != n_orbitals:
        raise ValueError(
            f"basis states span {basis_states.shape[1]} sites, expected "
            f"{n_orbitals}"
        )

    amplitudes = np.zeros(2 ** (2 * n_orbitals), dtype=complex)
    for row, coefficient in zip(basis_states, coefficients, strict=True):
        alpha, beta = alps_state_to_occupations(row)
        # Blocked order: modes 0..K-1 are alpha, K..2K-1 are beta. Little-endian
        # index: mode q contributes bit q.
        index = 0
        for orbital in range(n_orbitals):
            if alpha[orbital]:
                index |= 1 << orbital
            if beta[orbital]:
                index |= 1 << (n_orbitals + orbital)
        amplitudes[index] = coefficient * interleaved_to_blocked_parity(alpha, beta)
    return amplitudes
