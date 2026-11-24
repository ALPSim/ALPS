# SPDX-License-Identifier: MIT
"""The steady-state test on directed-loop SSE runs of the Heisenberg model.

Part 1 starts the spin-1/2 Heisenberg antiferromagnet on a 32 x 32 square
lattice at T = 0.25 without thermalization, several times. The energy relaxes
within the first bin: a large but short transient, which the linear test
barely notices. Part 2 repeats an equilibrated run many times and records how
often the test passes.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import pyalps


def steady_state_z(ts):
    """Return the fitted slope, its standard error and z = |slope| / error."""
    ts = np.asarray(ts, dtype=float)
    N = len(ts)
    t = np.arange(1, N + 1)
    slope = np.polyfit(t, ts, 1)[0]
    error = np.std(ts, ddof=1) * np.sqrt(12. / (N * (N**2 - 1)))
    return slope, error, abs(slope) / error


def cutoff(c):
    """z0 such that a standard normal Z satisfies P(|Z| < z0) = c."""
    return norm.ppf((1. + c) / 2.)


# Part 1: runs that start out of equilibrium. Without thermalization,
# dirloop_sse needs INITIAL_CUTOFF, the initial length of the operator string.
NSTART = 10
parms = [{
    'LATTICE'        : "square lattice",
    'MODEL'          : "spin",
    'local_S'        : 0.5,
    'L'              : 32,
    'J'              : 1.,
    'T'              : 0.25,
    'THERMALIZATION' : 0,
    'INITIAL_CUTOFF' : 40000,
    'SWEEPS'         : 2000,
} for run in range(NSTART)]
input_file = pyalps.writeInputFiles('parm1c_sse', parms)
pyalps.runApplication('dirloop_sse', input_file, Tmin=1)

files = pyalps.getResultFiles(prefix='parm1c_sse')
series = [np.array(pyalps.loadTimeSeries(f, 'Energy')) for f in files]
z = np.array([steady_state_z(ts)[2] for ts in series])
z_rest = np.array([steady_state_z(ts[1:])[2] for ts in series])
print('runs without thermalization: z =', ' '.join('%.1f' % v for v in z))
print('the same runs without the first bin: z =', ' '.join('%.1f' % v for v in z_rest))
for c in [0.6827, 0.95]:
    print('fraction passing at c = %g (z0 = %.2f): %.1f with the first bin, %.1f without'
          % (c, cutoff(c), np.mean(z < cutoff(c)), np.mean(z_rest < cutoff(c))))

# each time-series entry is the sum of the energy over one bin
ts = series[0]
E = pyalps.loadMeasurements(files[:1], 'Energy')[0][0].y[0].mean
binsize = 2 ** int(round(np.log2(ts.mean() / E)))
sweeps = (np.arange(len(ts)) + 0.5) * binsize
plt.figure()
plt.plot(sweeps, ts / binsize, 'o-', markersize=3)
plt.xlabel('Monte Carlo sweeps (no thermalization)')
plt.ylabel('Energy')
plt.title('32 x 32 Heisenberg antiferromagnet, T = 0.25: z = %.1f' % z[0])
plt.show()

# Part 2: many equilibrated runs
NRUNS = 100
c = 0.6827
parms = [{
    'LATTICE'        : "square lattice",
    'MODEL'          : "spin",
    'local_S'        : 0.5,
    'L'              : 32,
    'J'              : 1.,
    'T'              : 0.25,
    'THERMALIZATION' : 1000,
    'SWEEPS'         : 4000,
} for run in range(NRUNS)]
input_file = pyalps.writeInputFiles('parm1c_sse_equilibrated', parms)
pyalps.runApplication('dirloop_sse', input_file, Tmin=1)

z = np.array([steady_state_z(pyalps.loadTimeSeries(f, 'Energy'))[2]
              for f in pyalps.getResultFiles(prefix='parm1c_sse_equilibrated')])
print('%d equilibrated runs, fraction passing at c = %g: %.2f (expected %.2f)'
      % (len(z), c, np.mean(z < cutoff(c)), c))

plt.figure()
plt.hist(z, bins=np.linspace(0, 5, 26), density=True, histtype='step', lw=2, label='SSE energy')
x = np.linspace(0, 5, 201)
plt.plot(x, 2 * norm.pdf(x), 'k--', label='|Z|, Z standard normal')
plt.axvline(cutoff(c), color='gray', lw=1)
plt.xlabel('z = |slope| / error')
plt.ylabel('probability density')
plt.title('Equilibrated runs, 32 x 32 Heisenberg antiferromagnet, T = 0.25')
plt.legend()
plt.show()
