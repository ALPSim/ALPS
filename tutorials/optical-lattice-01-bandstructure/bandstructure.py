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
"""

import numpy as np

# Physical constants (SI)
h    = 6.62607015e-34
hbar = h / (2 * np.pi)
kB   = 1.380649e-23
amu  = 1.66053906660e-27
bohr = 5.29177210903e-11

V0   = np.array([8., 8., 8.])        # lattice depth in recoil energies
wlen = np.array([843., 843., 843.])  # laser wavelength in nanometer
a    = 114.8                         # s-wave scattering length in bohr radius
m    = 86.99                         # mass in atomic mass unit
L    = 200                           # lattice size (along 1 direction)
M    = 20                            # plane waves e^{i 2 m pi x}, m = -M..M


def trapezoid(f, dx):
    """Trapezoidal rule for samples f on a uniform grid with spacing dx."""
    return dx * (f.sum() - (f[0] + f[-1]) / 2)


def band_1d(V0, L, M):
    """Lowest band of one lattice direction.

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

    # w(x) decays exponentially, so a few lattice spacings around its site suffice.
    x, dx = np.linspace(-4, 4, 4001, retstep=True)
    q = ks[:, None] + ms[None, :]                          # q = k_x + m
    w = (c[:, :, None] * np.exp(2j * np.pi * q[:, :, None] * x)).sum((0, 1)).real / L
    return t, w, dx


Er2nK = h**2 / (2 * m * amu * (wlen * 1e-9)**2) / kB * 1e9

t, norm, w4 = np.empty(3), np.empty(3), np.empty(3)
for d in range(3):
    t[d], w, dx = band_1d(V0[d], L, M)
    norm[d] = trapezoid(w**2, dx)
    w4[d] = trapezoid(w**4, dx) / (wlen[d] / 2 * 1e-9)    # in 1/m

t_nK = t * Er2nK
U_nK = 4 * np.pi * hbar**2 * a * bohr / (m * amu) * np.prod(w4) / kB * 1e9

print(f'Er2nK  = {Er2nK}')
print(f't [nK] = {t_nK}')
print(f'U [nK] = {U_nK:.6g}')
print(f'U/t    = {U_nK / t_nK}')
print(f'norm of w (should be 1) = {norm}')
