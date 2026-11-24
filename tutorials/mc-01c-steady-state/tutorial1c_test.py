# SPDX-License-Identifier: MIT
"""The steady-state test of pyalps.checkSteadyState, computed by hand.

Defines the test statistic z and the cutoff z0, prints z0 for several
confidence levels, plots the acceptance region on the standard normal
density, and checks how often series without any drift pass: independent
samples pass at the nominal rate, autocorrelated samples less often.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


def steady_state_z(ts):
    """Return the fitted slope, its standard error and z = |slope| / error.

    This is the statistic of pyalps.checkSteadyState: a least-squares line
    through the points (1, ts[0]), ..., (N, ts[N-1]), with the slope error
    of N independent samples whose standard deviation is that of the series.
    """
    ts = np.asarray(ts, dtype=float)
    N = len(ts)
    t = np.arange(1, N + 1)
    slope = np.polyfit(t, ts, 1)[0]
    error = np.std(ts, ddof=1) * np.sqrt(12. / (N * (N**2 - 1)))
    return slope, error, abs(slope) / error


def cutoff(c):
    """z0 such that a standard normal Z satisfies P(|Z| < z0) = c."""
    return norm.ppf((1. + c) / 2.)


def cutoff_pyalps_300(c):
    """The cutoff of pyalps 3.0.0 and earlier, which used the level 1 - c."""
    return norm.ppf(1. - c / 2.)


# The cutoff for common confidence levels, and how often a series without
# drift passes with each cutoff
print('     c      z0   pass rate   z0 (3.0.0)   pass rate (3.0.0)')
for c in [0.6827, 0.9, 0.95, 0.99]:
    z0, z0_old = cutoff(c), cutoff_pyalps_300(c)
    print('%6.4f   %5.3f   %9.3f   %10.3f   %17.3f'
          % (c, z0, 2 * norm.cdf(z0) - 1, z0_old, 2 * norm.cdf(z0_old) - 1))

# The acceptance region |z| < z0 on the standard normal density
x = np.linspace(-4, 4, 801)
plt.figure()
plt.plot(x, norm.pdf(x), color='black')
for c, color in [(0.95, 'tab:blue'), (0.6827, 'tab:orange')]:
    z0 = cutoff(c)
    inside = np.abs(x) < z0
    plt.fill_between(x[inside], norm.pdf(x[inside]), color=color, alpha=0.35,
                     label='c = %g: |z| < %.2f' % (c, z0))
z0_old = cutoff_pyalps_300(0.6827)
plt.axvline(-z0_old, color='tab:red', ls='--')
plt.axvline(z0_old, color='tab:red', ls='--', label='pyalps 3.0.0, c = 0.6827: |z| < %.2f' % z0_old)
plt.xlabel('z = slope / error')
plt.ylabel('standard normal density')
plt.title('Acceptance region of the steady-state test')
plt.legend(loc='upper left', fontsize='small')
plt.show()

# Series without drift: independent Gaussian noise, and autocorrelated noise
# x[i] = a x[i-1] + noise with integrated autocorrelation time
# tau = (1 + a) / (2 (1 - a)), as in MC-01a
rng = np.random.default_rng(2026)
N, samples = 128, 10000
z0 = cutoff(0.6827)
print()
print('    tau   fraction passing at c = 0.6827 (N = %d points)' % N)
for a in [0., 0.5, 0.8, 0.95]:
    z = np.empty(samples)
    for k in range(samples):
        noise = rng.normal(size=N)
        series = np.empty(N)
        series[0] = noise[0] / np.sqrt(1 - a**2)
        for i in range(1, N):
            series[i] = a * series[i - 1] + noise[i]
        z[k] = steady_state_z(series)[2]
    tau = (1 + a) / (2 * (1 - a))
    print('%7.1f   %.3f' % (tau, np.mean(z < z0)))
