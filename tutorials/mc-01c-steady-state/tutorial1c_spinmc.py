# SPDX-License-Identifier: MIT
"""The steady-state test on classical Monte Carlo runs of the 2D Ising model.

Part 1 runs a large lattice at T_c without thermalization, so that the
magnetization is still relaxing from the ordered initial state: the test
detects this drift. Part 2 repeats an equilibrated run many times with local
and with cluster updates and records how often the test passes, which shows
how autocorrelation between the time-series bins inflates the false-alarm
rate.
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


# Part 1: a run that is still drifting
parms = [{
    'LATTICE'        : "square lattice",
    'MODEL'          : "Ising",
    'L'              : 128,
    'J'              : 1.,
    'T'              : 2.269186,
    'UPDATE'         : "local",
    'THERMALIZATION' : 0,
    'SWEEPS'         : 2000,
}]
input_file = pyalps.writeInputFiles('parm1c_drift', parms)
pyalps.runApplication('spinmc', input_file, Tmin=1, writexml=True)

files = pyalps.getResultFiles(prefix='parm1c_drift')
ts = np.array(pyalps.loadTimeSeries(files[0], '|Magnetization|'))
slope, error, z = steady_state_z(ts)
print('drift run: %d bins, z = %.2f, cutoff z0 = %.2f, passes: %s'
      % (len(ts), z, cutoff(0.6827), z < cutoff(0.6827)))

data = pyalps.loadMeasurements(files, '|Magnetization|')
for d in pyalps.checkSteadyState(data):
    s = d.props['checkSteadyState']
    print('pyalps.checkSteadyState: z = %.2f, z0 = %.2f, passes: %s'
          % (s['statistics']['z'], s['statistics']['z0'], s['value']))

# each time-series entry is the sum of |M| per site over one bin
M = data[0][0].y[0].mean
binsize = 2 ** int(round(np.log2(ts.mean() / M)))
sweeps = (np.arange(len(ts)) + 0.5) * binsize
plt.figure()
plt.plot(sweeps, ts / binsize, 'o-', markersize=3)
plt.xlabel('Monte Carlo sweeps (no thermalization)')
plt.ylabel('|Magnetization| per site')
plt.title('L = 128, T = T_c, local updates: z = %.1f' % z)
plt.show()

# Part 2: many equilibrated runs, with local and with cluster updates
NRUNS = 100
c = 0.6827
plt.figure()
for update in ['local', 'cluster']:
    parms = [{
        'LATTICE'        : "square lattice",
        'MODEL'          : "Ising",
        'L'              : 48,
        'J'              : 1.,
        'T'              : 2.269186,
        'UPDATE'         : update,
        'THERMALIZATION' : 10000,
        'SWEEPS'         : 20000,
    } for run in range(NRUNS)]
    prefix = 'parm1c_' + update
    input_file = pyalps.writeInputFiles(prefix, parms)
    pyalps.runApplication('spinmc', input_file, Tmin=1, writexml=True)

    z = np.array([steady_state_z(pyalps.loadTimeSeries(f, '|Magnetization|'))[2]
                  for f in pyalps.getResultFiles(prefix=prefix)])
    print('%-7s updates: %d runs, fraction passing at c = %g: %.2f (expected %.2f)'
          % (update, len(z), c, np.mean(z < cutoff(c)), c))
    plt.hist(z, bins=np.linspace(0, 5, 26), density=True, histtype='step', lw=2,
             label='%s updates' % update)

x = np.linspace(0, 5, 201)
plt.plot(x, 2 * norm.pdf(x), 'k--', label='|Z|, Z standard normal')
plt.axvline(cutoff(c), color='gray', lw=1)
plt.xlabel('z = |slope| / error')
plt.ylabel('probability density')
plt.title('Equilibrated runs, L = 48, T = T_c')
plt.legend()
plt.show()
