"""Dense Jordan-Wigner reference: eval ALPS-style term strings as matrices."""
import numpy as np


class Op:
    def __init__(self, m): self.m = m
    def __add__(self, o): return Op(self.m + (o.m if isinstance(o, Op) else o * np.eye(len(self.m))))
    __radd__ = __add__
    def __sub__(self, o): return self + (-1) * o
    def __rsub__(self, o): return (-1) * self + o
    def __neg__(self): return (-1) * self
    def __mul__(self, o): return Op(self.m @ o.m) if isinstance(o, Op) else Op(self.m * o)
    __rmul__ = lambda self, o: Op(o * self.m)


def fock_ops(nsites):
    M = 2 * nsites
    D = 1 << M
    def ann(mode):
        A = np.zeros((D, D))
        for s in range(D):
            if s >> mode & 1:
                A[s ^ (1 << mode), s] = (-1) ** bin(s & ((1 << mode) - 1)).count("1")
        return A
    ops = {}
    for spin, a in (("up", 0), ("down", 1)):
        c = [ann(2 * i + a) for i in range(nsites)]
        ops[f"c_{spin}"] = lambda i, c=c: Op(c[i])
        ops[f"cdag_{spin}"] = lambda i, c=c: Op(c[i].T.copy())
    return ops


def dense(nsites, bonds, plaquettes, bond_expr, plaq_expr, roles, params):
    ops = fock_ops(nsites)
    H = Op(np.zeros((1 << 2 * nsites,) * 2))
    for i, j in bonds:
        H = H + eval(bond_expr, {}, {**ops, **params, "i": i, "j": j})
    for sites in plaquettes:
        H = H + eval(plaq_expr, {}, {**ops, **params, **dict(zip(roles, sites))})
    assert np.allclose(H.m, H.m.T)
    return np.linalg.eigvalsh(H.m)


