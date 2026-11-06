"""Dense Jordan-Wigner reference: eval ALPS-style term strings as matrices."""
import numpy as np
import scipy.sparse as sp


class Op:
    def __init__(self, m): self.m = m
    def __add__(self, o):
        if isinstance(o, Op):
            return Op(self.m + o.m)
        eye = sp.identity(self.m.shape[0], format="csr") if sp.issparse(self.m) else np.eye(len(self.m))
        return Op(self.m + o * eye)
    __radd__ = __add__
    def __sub__(self, o): return self + (-1) * o
    def __rsub__(self, o): return (-1) * self + o
    def __neg__(self): return (-1) * self
    def __mul__(self, o): return Op(self.m @ o.m) if isinstance(o, Op) else Op(self.m * o)
    __rmul__ = lambda self, o: Op(o * self.m)


def fock_ops(nsites, sparse=False):
    """c and cdag of every site as Op; scipy sparse matrices if sparse."""
    M = 2 * nsites
    D = 1 << M
    def ann(mode):
        rows, cols, vals = [], [], []
        for s in range(D):
            if s >> mode & 1:
                rows.append(s ^ (1 << mode))
                cols.append(s)
                vals.append((-1) ** bin(s & ((1 << mode) - 1)).count("1"))
        A = sp.csr_matrix((vals, (rows, cols)), shape=(D, D))
        return A if sparse else A.toarray()
    ops = {}
    for spin, a in (("up", 0), ("down", 1)):
        c = [ann(2 * i + a) for i in range(nsites)]
        ops[f"c_{spin}"] = lambda i, c=c: Op(c[i])
        ops[f"cdag_{spin}"] = lambda i, c=c: Op(c[i].T.copy())
    return ops


def dense(nsites, bonds, groups, bond_expr, group_expr, roles, params):
    ops = fock_ops(nsites)
    H = Op(np.zeros((1 << 2 * nsites,) * 2))
    for i, j in bonds:
        H = H + eval(bond_expr, {}, {**ops, **params, "i": i, "j": j})
    for sites in groups:
        H = H + eval(group_expr, {}, {**ops, **params, **dict(zip(roles, sites))})
    assert np.allclose(H.m, H.m.T)
    return np.linalg.eigvalsh(H.m)


