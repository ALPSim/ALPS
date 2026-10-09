# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
"""Boson Hubbard parameters t and U for a cubic optical lattice.

Solves the lowest band of the separable lattice potential
V(r) = sum_alpha V0_alpha sin^2(pi x_alpha) in the plane-wave basis, builds the
Wannier function in each direction and evaluates

    t_alpha = -(1/L) sum_k eps_k cos(2 pi k)
    U       = g prod_alpha int |w(x_alpha)|^4 dx_alpha,   g = 4 pi hbar^2 a / m

Energies are in recoil energies E_r = h^2 / (2 m lambda^2) and lengths in lattice
spacings lambda/2 until converted to nK at the end. This replaces
pyalps.dwa.bandstructure, which was removed together with the DWA application.

Usage, from this directory:

    import bandstructure
    t, U = bandstructure.hubbard_parameters(V0, wlen, a, m, L)
"""

import numpy as np

# Physical constants (SI)
h    = 6.62607015e-34
hbar = h / (2 * np.pi)
kB   = 1.380649e-23
amu  = 1.66053906660e-27
bohr = 5.29177210903e-11


def trapezoid(f, dx):
    """Trapezoidal rule for samples f on a uniform grid with spacing dx."""
    return dx * (f.sum() - (f[0] + f[-1]) / 2)


def band_1d(V0, L, M):
    """Lowest band of one lattice direction of depth V0 (in E_r).

    Returns the hopping t in E_r, and the Wannier function w sampled on a grid x
    in lattice spacings with spacing dx.
    """
    ks = (np.arange(L) - L // 2) / L                       # k_x in [-1/2, 1/2)
    ms = np.arange(-M, M + 1)
    eps = np.empty(L)
    c = np.empty((L, ms.size))
    for i, k in enumerate(ks):
        H = np.diag(4. * (ms + k)**2 + V0 / 2.) \
            - V0 / 4. * (np.eye(ms.size, k=1) + np.eye(ms.size, k=-1))
        e, v = np.linalg.eigh(H)
        eps[i], c[i] = e[0], v[:, 0] * np.sign(v[M, 0])   # fix the gauge: c_0 > 0
    t = -np.mean(eps * np.cos(2 * np.pi * ks))

    # w(x) decays exponentially, so ten lattice spacings around its site suffice.
    # With c_0 > 0 the Wannier function is real and even, a sum of cosines.
    x, dx = np.linspace(-10, 10, 10001, retstep=True)
    w = np.zeros_like(x)
    for i, k in enumerate(ks):
        w += c[i] @ np.cos(2 * np.pi * np.outer(ms + k, x))
    return t, w / L, dx


def hubbard_parameters(V0, wlen, a, m, L, M=20):
    """Hopping t (one per direction) and onsite interaction U, both in nK.

    V0   -- lattice depth in recoil energies, one value per direction
    wlen -- laser wavelength in nanometer, one value per direction
    a    -- s-wave scattering length in Bohr radii
    m    -- mass in atomic mass units
    L    -- lattice size along one direction
    M    -- plane-wave cutoff, e^{i 2 m pi x} with m = -M..M
    """
    V0, wlen = np.asarray(V0, dtype=float), np.asarray(wlen, dtype=float)
    Er2nK = h**2 / (2 * m * amu * (wlen * 1e-9)**2) / kB * 1e9

    t, w4 = np.empty(len(V0)), np.empty(len(V0))
    for d in range(len(V0)):
        t[d], w, dx = band_1d(V0[d], L, M)
        w4[d] = trapezoid(w**4, dx) / (wlen[d] / 2 * 1e-9)    # in 1/m

    U = 4 * np.pi * hbar**2 * a * bohr / (m * amu) * np.prod(w4) / kB * 1e9
    return t * Er2nK, U
