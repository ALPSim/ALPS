# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
"""Density profile of trapped bosons in a cubic optical lattice.

Simulates the Bose-Hubbard model with a parabolic trap,

    H = -t sum_<ij> b_i^+ b_j + U/2 sum_i n_i (n_i - 1) - sum_i (mu - V_T(r_i)) n_i,

with the directed-loop SSE quantum Monte Carlo code and plots the local density
<n_i>. U/t = 8.11 is the value optical-lattice-01-bandstructure gives for a
V0 = 8 E_r lattice of 87Rb. This replaces the dwa-02-density-profile tutorial,
which was removed together with the DWA application.

The expected output is shown in density_profile.png. The lattice is kept small
so the tutorial runs in a few minutes, which means it still sits inside the
cloud: the density falls off toward the edges but does not quite vanish
there. A bigger lattice (larger L) lets the density drop closer to zero at
the edges.
"""

import numpy as np
import matplotlib.pyplot as plt
import pyalps

L = 9
K = 0.65                                        # trap curvature V_T = K r^2, in units of t
c = (L - 1) / 2.

parms = [{
    'LATTICE' : 'inhomogeneous simple cubic lattice',
    'L'       : L,
    'MODEL'   : 'boson Hubbard',
    'Nmax'    : 4,
    't'       : 1.,
    'U'       : 8.11,
    'mu'      : '4.05 - %g*((x-%g)*(x-%g) + (y-%g)*(y-%g) + (z-%g)*(z-%g))' % ((K,) + (c,) * 6),
    'T'       : 1.,
    'THERMALIZATION' : 1000,
    'SWEEPS'         : 5000,
    'MEASURE_LOCAL[Local Density]' : 'n',
}]

input_file = pyalps.writeInputFiles('parm_trap', parms)
pyalps.runApplication('dirloop_sse', input_file)

data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm_trap'), 'Local Density')[0][0]
n = np.asarray(data.y.mean).reshape(L, L, L)
print('Total number of bosons: %.2f' % n.sum())
print('Density at the trap center: %.3f' % n[L // 2, L // 2, L // 2])

r = np.arange(L) - c
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(r, n[:, L // 2, L // 2], 'o-')
ax1.set_xlabel('$x$ (lattice sites from trap center)')
ax1.set_ylabel(r'$\langle n_i \rangle$')
ax1.set_title('Cut through the trap center')
im = ax2.imshow(n[:, :, L // 2], extent=[r[0], r[-1], r[0], r[-1]], origin='lower')
fig.colorbar(im, ax=ax2, label=r'$\langle n_i \rangle$')
ax2.set_title('Center layer $z = 0$')
plt.tight_layout()
plt.show()
