"""Exact reference states in ALPS's own local basis, to compare with dmrg's MPS.

The MPS from dmrg SAVE_MPS=1 (ALPSim/ALPS#216) is indexed by the ALPS site
basis (order and quantum numbers in /mps/site_basis). The reference here is
built in that same basis: site operators as matrices on each local state,
embedded with Kronecker products (site 0 most significant, as the MPS is
contracted). Fermion operators carry the Jordan-Wigner string over the sites
before theirs, and on one site up comes before down, which is the ALPS
convention for products of fermion operators.
"""
import numpy as np
import scipy.sparse as sparse


class Op:
    """Sparse operator with the arithmetic of ALPS expressions."""
    def __init__(self, m): self.m = m
    def __add__(self, o):
        if isinstance(o, Op):
            return Op(self.m + o.m)
        return Op(self.m + o * sparse.identity(self.m.shape[0], format="csr"))
    __radd__ = __add__
    def __sub__(self, o): return self + (-1) * o
    def __rsub__(self, o): return (-1) * self + o
    def __neg__(self): return (-1) * self
    def __mul__(self, o): return Op(self.m @ o.m) if isinstance(o, Op) else Op(self.m * o)
    def __rmul__(self, o): return Op(o * self.m)


def local_operators(basis, S=None):
    """Single-site matrices [out, in] for the quantum numbers of each local state."""
    d = len(next(iter(basis.values())))
    ops = {}

    def mat(f):
        m = np.zeros((d, d))
        for i in range(d):
            for j in range(d):
                m[j, i] = f(i, j)
        return m

    if "Nup" in basis:                                  # spinful fermions
        nu, nd = basis["Nup"].astype(int), basis["Ndown"].astype(int)
        ops["c_up"] = mat(lambda i, j: float(nu[i] == 1 and nu[j] == 0 and nd[i] == nd[j]))
        ops["c_down"] = mat(lambda i, j: (-1.) ** nu[i] * (nd[i] == 1 and nd[j] == 0 and nu[i] == nu[j]))
        ops["cdag_up"], ops["cdag_down"] = ops["c_up"].T, ops["c_down"].T
        ops["n_up"], ops["n_down"] = np.diag(nu * 1.), np.diag(nd * 1.)
        parity = np.diag((-1.) ** (nu + nd))
        return ops, parity, ("c_up", "c_down", "cdag_up", "cdag_down")
    if "Sz" in basis:                                   # spins
        m = basis["Sz"].astype(float)
        s = float(basis["S"][0]) if "S" in basis else (d - 1) / 2 if S is None else S
        ops["Sz"] = np.diag(m)
        ops["Splus"] = mat(lambda i, j: np.sqrt(s * (s + 1) - m[i] * (m[i] + 1)) * np.isclose(m[j], m[i] + 1))
        ops["Sminus"] = ops["Splus"].T
        return ops, None, ()
    if "N" in basis:                                    # bosons
        n = basis["N"].astype(int)
        ops["b"] = mat(lambda i, j: np.sqrt(n[i]) * (n[j] == n[i] - 1))
        ops["bdag"] = ops["b"].T
        ops["n"] = np.diag(n * 1.)
        return ops, None, ()
    raise ValueError(f"unknown site basis {sorted(basis)}")


def site_operators(nsites, basis, S=None):
    """{name: function(site) -> Op on the full space}."""
    loc, parity, fermionic = local_operators(basis, S)
    d = len(next(iter(basis.values())))
    eye = sparse.identity(d, format="csr")

    def embed(name, i):
        m = sparse.identity(1, format="csr")
        for k in range(nsites):
            if k == i:
                f = sparse.csr_matrix(loc[name])
            elif k < i and name in fermionic:
                f = sparse.csr_matrix(parity)
            else:
                f = eye
            m = sparse.kron(m, f, format="csr")
        return Op(m)

    return {name: (lambda i, name=name: embed(name, i)) for name in loc}


def totals(nsites, basis):
    """Total of each quantum number for every product state."""
    out = {}
    for q, v in basis.items():
        t = np.zeros(1)
        for _ in range(nsites):
            t = (t[:, None] + v[None, :]).ravel()
        out[q] = t
    return out


def hamiltonian(nsites, basis, terms, params=None, S=None):
    """terms: list of (expression, [site tuples]); roles bound in order."""
    ops = site_operators(nsites, basis, S)
    d = len(next(iter(basis.values())))
    H = Op(sparse.csr_matrix((d ** nsites,) * 2))
    for expr, roles, groups in terms:
        for g in groups:
            H = H + eval(expr, {}, {**ops, **(params or {}), **dict(zip(roles, g))})
    assert abs(H.m - H.m.T).max() < 1e-10
    return H.m


def ground_state(H, nsites, basis, sector):
    """Lowest eigenpair in a sector {quantum number: total}; the full state vector."""
    tot = totals(nsites, basis)
    keep = np.ones(H.shape[0], bool)
    for q, v in sector.items():
        keep &= np.isclose(tot[q], v)
    idx = np.nonzero(keep)[0]
    e, v = np.linalg.eigh(H[idx][:, idx].toarray())
    psi = np.zeros(H.shape[0])
    psi[idx] = v[:, 0]
    return e, psi


def contract(tensors):
    """Full state vector of an MPS, site 0 most significant."""
    psi = tensors[0][0]                                  # [s0, r]
    for t in tensors[1:]:
        psi = np.tensordot(psi, t, axes=([-1], [0]))
    return psi[..., 0].ravel()


def compare(mps, ref_psi, H):
    """(|overlap|, <H> of the MPS, max deviation of the components after fixing the sign)."""
    psi = contract(mps["tensors"])
    norm = np.linalg.norm(psi)
    psi = psi / norm
    ov = float(np.dot(ref_psi, psi))
    energy = float(psi @ (H @ psi))
    dev = float(np.max(np.abs(np.sign(ov) * psi - ref_psi)))
    return abs(ov), energy, dev, norm
