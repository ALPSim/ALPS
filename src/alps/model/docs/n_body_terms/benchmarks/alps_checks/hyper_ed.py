"""Exact diagonalization of spinful fermion models with n-site terms.

Pipeline:
    TOML model (term strings on hyperedges)
      -> parse + normal-order       -> coefficient tables h_ab, V_abcd, ...
      -> fill Fock-space H          -> sparse D x D matrix in an (Nup, Ndown) sector
      -> diagonalize                -> energies, ground state, observables

A term is attached to a class of hyperedges ("site", "bond", "plaquette", ...),
each an ordered tuple of sites, so 1-, 2- and 4-site terms are handled by the
same code. Modes are numbered 2*site + spin (up = 0, down = 1); Jordan-Wigner
signs follow that order.

Usage:  python hyper_ed.py model.toml
"""

import ast
from collections import defaultdict
from itertools import combinations
import sys
import tomllib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

UP, DOWN = 0, 1


# ---------------------------------------------------------------- operator algebra

class OpPoly:
    """A polynomial in fermion operators: {(op, op, ...): coeff}, op = (dagger, mode)."""

    def __init__(self, terms=None):
        self.terms = dict(terms or {})

    @staticmethod
    def scalar(x):
        return OpPoly({(): x})

    def _lift(self, other):
        return other if isinstance(other, OpPoly) else OpPoly.scalar(other)

    def __add__(self, other):
        out = defaultdict(complex, self.terms)
        for word, c in self._lift(other).terms.items():
            out[word] += c
        return OpPoly(out)

    __radd__ = __add__

    def __neg__(self):
        return OpPoly({w: -c for w, c in self.terms.items()})

    def __sub__(self, other):
        return self + (-self._lift(other))

    def __rsub__(self, other):
        return self._lift(other) - self

    def __mul__(self, other):
        other = self._lift(other)
        out = defaultdict(complex)
        for w1, c1 in self.terms.items():
            for w2, c2 in other.terms.items():
                out[w1 + w2] += c1 * c2
        return OpPoly(out)

    def __rmul__(self, other):
        return self._lift(other) * self

    def __truediv__(self, x):
        return self * (1.0 / x)


def as_poly(x):
    return x if isinstance(x, OpPoly) else OpPoly.scalar(x)


def normal_order(word):
    """Return {(creators, annihilators): coeff} for an operator word.

    Uses {c_a, c+_b} = delta_ab to move creators left, then sorts each group
    by mode with the permutation sign. Repeated modes in a group vanish.
    """
    out = defaultdict(complex)
    stack = [(tuple(word), 1.0)]
    while stack:
        w, c = stack.pop()
        for k in range(len(w) - 1):
            (d1, m1), (d2, m2) = w[k], w[k + 1]
            if not d1 and d2:                      # c_m1 c+_m2 -> delta - c+_m2 c_m1
                if m1 == m2:
                    stack.append((w[:k] + w[k + 2:], c))
                stack.append((w[:k] + (w[k + 1], w[k]) + w[k + 2:], -c))
                break
        else:
            cre = [m for d, m in w if d]
            ann = [m for d, m in w if not d]
            if len(set(cre)) < len(cre) or len(set(ann)) < len(ann):
                continue
            sign = _sort_sign(cre) * _sort_sign(ann)
            out[(tuple(sorted(cre)), tuple(sorted(ann)))] += sign * c
    return out


def _sort_sign(modes):
    sign, modes = 1, list(modes)
    for i in range(len(modes)):
        for j in range(len(modes) - 1 - i):
            if modes[j] > modes[j + 1]:
                modes[j], modes[j + 1] = modes[j + 1], modes[j]
                sign = -sign
    return sign


# ---------------------------------------------------------------- term strings

def _site_operators(site):
    def op(dag, spin):
        return OpPoly({((dag, 2 * site + spin),): 1.0})
    cu, cd = op(False, UP), op(False, DOWN)
    cdu, cdd = op(True, UP), op(True, DOWN)
    nu, nd = cdu * cu, cdd * cd
    return {
        "c_up": cu, "c_down": cd, "cdag_up": cdu, "cdag_down": cdd,
        "n_up": nu, "n_down": nd, "n": nu + nd,
        "Sz": 0.5 * (nu - nd), "Splus": cdu * cd, "Sminus": cdd * cu,
    }


def evaluate(expr, roles, params, site_ops):
    """Evaluate a term string. roles maps role names (i, j, ...) to lattice sites;
    site_ops(site) returns {operator name: object} in the chosen algebra."""
    tree = ast.parse(expr, mode="eval").body

    def ev(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return node.value
        if isinstance(node, ast.Name):
            if node.id not in params:
                raise ValueError(f"unknown parameter {node.id!r} in {expr!r}")
            return params[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            v = ev(node.operand)
            return -v if isinstance(node.op, ast.USub) else v
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
            a, b = ev(node.left), ev(node.right)
            return {ast.Add: lambda: a + b, ast.Sub: lambda: a - b,
                    ast.Mult: lambda: a * b, ast.Div: lambda: a / b}[type(node.op)]()
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and len(node.args) == 1 and isinstance(node.args[0], ast.Name)):
            role = node.args[0].id
            if role not in roles:
                raise ValueError(f"site {role!r} is not one of {list(roles)} in {expr!r}")
            ops = site_ops(roles[role])
            if node.func.id not in ops:
                raise ValueError(f"unknown operator {node.func.id!r} in {expr!r}")
            return ops[node.func.id]
        raise ValueError(f"unsupported syntax in {expr!r}: {ast.dump(node)}")

    return ev(tree)


# ---------------------------------------------------------------- lattices as hypergraphs

def hypergraph(lattice):
    """Return (number of sites, {edge class: [ordered site tuples]})."""
    kind, L = lattice["kind"], lattice.get("L", 0)
    periodic = lattice.get("boundary", "open") == "periodic"
    edges = defaultdict(list)
    if kind == "chain":
        n = L
        edges["site"] = [(i,) for i in range(L)]
        edges["bond"] = [(i, i + 1) for i in range(L - 1)]
        if periodic and L > 2:
            edges["bond"].append((L - 1, 0))
    elif kind == "ladder":                     # site (x, y) -> 2x + y, y = 0, 1
        n = 2 * L
        s = lambda x, y: 2 * (x % L) + y
        xs = range(L) if periodic and L > 2 else range(L - 1)
        edges["site"] = [(i,) for i in range(n)]
        edges["rung"] = [(s(x, 0), s(x, 1)) for x in range(L)]
        edges["leg"] = [(s(x, y), s(x + 1, y)) for x in xs for y in (0, 1)]
        edges["bond"] = edges["rung"] + edges["leg"]
        # plaquette roles (i, j, k, l): rung x = (i, j), rung x+1 = (l, k)
        edges["plaquette"] = [(s(x, 0), s(x, 1), s(x + 1, 1), s(x + 1, 0)) for x in xs]
    elif kind == "custom":
        n = lattice["sites"]
    else:
        raise ValueError(f"unknown lattice kind {kind!r}")
    for cls, tuples in lattice.get("edges", {}).items():
        edges[cls] += [tuple(t) for t in tuples]
    return n, dict(edges)


# ---------------------------------------------------------------- model -> coefficients

def coefficients(terms, n_sites, edges, params):
    """Sum every term over its hyperedges and normal-order the result.

    Returns {(creators, annihilators): coeff}; entries with one creator and one
    annihilator are h_ab, two and two are V_abcd, and so on.
    """
    table = defaultdict(complex)
    for term in terms:
        roles = term["sites"]
        cls = term["on"]
        if cls not in edges:
            raise ValueError(f"lattice has no {cls!r} hyperedges")
        for edge in edges[cls]:
            if len(edge) != len(roles):
                raise ValueError(f"{cls!r} edges have {len(edge)} sites, term names {roles}")
            poly = evaluate(term["expr"], dict(zip(roles, edge)), params, _site_operators)
            for word, c in as_poly(poly).terms.items():
                for key, c2 in normal_order(word).items():
                    table[key] += c * c2
    return {k: v for k, v in table.items() if abs(v) > 1e-12}


# ---------------------------------------------------------------- Fock space

def sector_basis(n_sites, n_up, n_down):
    """Bit-string basis states with fixed Nup and Ndown, and their index.

    Lin ordering: the up arrangements and the down arrangements are each
    numbered 0, 1, 2, ... in the order itertools.combinations lists them, and
    index = (up number) * (number of down arrangements) + (down number).
    The bit string (bit 2*site + spin set if occupied) is kept for the signs.
    """
    downs_list = list(combinations(range(n_sites), n_down))
    states = []
    for ups in combinations(range(n_sites), n_up):
        up_bits = sum(1 << (2 * s + UP) for s in ups)
        for downs in downs_list:
            states.append(up_bits | sum(1 << (2 * s + DOWN) for s in downs))
    return states, {s: i for i, s in enumerate(states)}


def apply_word(state, cre, ann):
    """Apply c+_cre[0] ... c+_cre[-1] c_ann[0] ... c_ann[-1] to a basis state.
    Returns (new state, sign) or None if the result is zero."""
    sign = 1
    for m in reversed(ann):                   # rightmost operator acts first
        if not state >> m & 1:
            return None
        if bin(state & ((1 << m) - 1)).count("1") & 1:
            sign = -sign
        state ^= 1 << m
    for m in reversed(cre):
        if state >> m & 1:
            return None
        if bin(state & ((1 << m) - 1)).count("1") & 1:
            sign = -sign
        state ^= 1 << m
    return state, sign


def fock_matrix(table, states, index):
    rows, cols, vals = [], [], []
    for col, s in enumerate(states):
        for (cre, ann), c in table.items():
            res = apply_word(s, cre, ann)
            if res is None:
                continue
            t, sign = res
            if t not in index:
                raise ValueError("term leaves the (Nup, Ndown) sector")
            rows.append(index[t])
            cols.append(col)
            vals.append(sign * c)
    D = len(states)
    H = sp.csr_matrix((vals, (rows, cols)), shape=(D, D), dtype=complex)
    if np.allclose(H.imag.data, 0) if H.nnz else True:
        H = H.real
    return H


def lowest(H, k):
    D = H.shape[0]
    if D <= 2000:
        w, v = np.linalg.eigh(H.toarray())
        return w[:k], v[:, :k]
    w, v = sla.eigsh(H, k=k, which="SA")
    order = np.argsort(w)
    return w[order], v[:, order]


# ---------------------------------------------------------------- driver

def load(path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def solve(model, params=None, n_levels=4, verbose=True):
    params = {**model.get("parameters", {}), **(params or {})}
    n_sites, edges = hypergraph(model["lattice"])
    table = coefficients(model["term"], n_sites, edges, params)
    states, index = sector_basis(n_sites, model["sector"]["Nup"], model["sector"]["Ndown"])
    H = fock_matrix(table, states, index)
    if abs(H - H.conj().T).max() > 1e-10:
        raise ValueError("Hamiltonian is not Hermitian; check the term strings")
    energies, vectors = lowest(H, min(n_levels, len(states)))
    if verbose:
        by_rank = defaultdict(int)
        for cre, ann in table:
            by_rank[f"{len(cre)}+{len(ann)}"] += 1
        print(f"sites: {n_sites}   modes: {2 * n_sites}   hyperedges: "
              + ", ".join(f"{k}={len(v)}" for k, v in edges.items()))
        print("coefficients (creators+annihilators: count): "
              + ", ".join(f"{k}: {v}" for k, v in sorted(by_rank.items())))
        print(f"sector Nup={model['sector']['Nup']} Ndown={model['sector']['Ndown']}: "
              f"D = {len(states)}, nonzeros in H = {H.nnz}")
        print("lowest energies:", "  ".join(f"{e:.10f}" for e in energies))
    return {"energies": energies, "vectors": vectors, "H": H, "table": table,
            "states": states, "index": index, "n_sites": n_sites, "edges": edges,
            "params": params}


def expectation(result, expr, roles):
    """<psi0| O |psi0> for an operator string on the given sites."""
    poly = evaluate(expr, roles, result["params"], _site_operators)
    table = defaultdict(complex)
    for word, c in as_poly(poly).terms.items():
        for key, c2 in normal_order(word).items():
            table[key] += c * c2
    O = fock_matrix(table, result["states"], result["index"])
    psi = result["vectors"][:, 0]
    return (psi.conj() @ (O @ psi)).real


if __name__ == "__main__":
    solve(load(sys.argv[1]))
