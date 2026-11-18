#############################################################################
#
# ALPS Project Applications: Directed Worm Algorithm
#
# Copyright (C) 2013 by Lode Pollet      <pollet@phys.ethz.ch>
#                       Ping Nang Ma     <pingnang@phys.ethz.ch>
#                       Matthias Troyer  <troyer@phys.ethz.ch>
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
#
#############################################################################

# The headers
import pyalps
import numpy as np
import matplotlib.pyplot as plt

# Set up a python list of parameters (python) dictionaries. UPDATE selects
# local single-spin updates; without it spinmc uses cluster updates:
parms = [{
  'LATTICE'         : "square lattice",
  'MODEL'           : "Ising",
  'L'               : 48,
  'J'               : 1.,
  'T'               : 2.269186,
  'THERMALIZATION'  : 10000,
  'SWEEPS'          : 50000,
  'UPDATE'          : "local",
}]

# and an otherwise identical set with cluster updates (parm1b) for comparison:
parms_cluster = [dict(parms[0], UPDATE="cluster")]

# Write into XML input files:
input_file = pyalps.writeInputFiles('parm1a',parms)
input_file_cluster = pyalps.writeInputFiles('parm1b',parms_cluster)

# and run the application spinmc. Tmin is the minimum time in seconds between
# checks whether a simulation has finished; the simulation keeps measuring
# until the next check, so a run can contain more than SWEEPS measurements:
pyalps.runApplication('spinmc', input_file, Tmin=10, writexml=True)
pyalps.runApplication('spinmc', input_file_cluster, Tmin=10, writexml=True)

for prefix, update in [('parm1a', 'local'), ('parm1b', 'cluster')]:
    # We first get the list of all result files via:
    files = pyalps.getResultFiles(prefix=prefix)

    # and then extract, say the timeseries of the |Magnetization| measurements.
    # The timeseries is stored as at most ~128 bins; each entry is the sum of
    # the measurements in one bin, and the bin size (a power of two) doubles
    # as the run grows:
    ts_M = np.array(pyalps.loadTimeSeries(files[0], '|Magnetization|'))
    M = pyalps.loadMeasurements(files, '|Magnetization|')[0][0].y[0].mean
    binsize = 2 ** int(round(np.log2(ts_M.mean() / M)))
    sweeps = (np.arange(len(ts_M)) + 0.5) * binsize

    # We can then visualize graphically |M| per site against the sweep number:
    plt.figure()
    plt.plot(sweeps, ts_M / binsize, 'o-')
    plt.axhline(M, ls='--', color='gray', label='mean')
    plt.ylim(0, 1)
    plt.xlabel('Monte Carlo sweep')
    plt.ylabel('|Magnetization| per site')
    plt.title('T = T_c, %s updates (each point: mean over %d sweeps)' % (update, binsize))
    plt.legend()
plt.show()

# With local updates |M| still swings strongly between bins: at T_c the
# autocorrelation time of single-spin updates is hundreds of sweeps or more
# (critical slowing down). With cluster updates the bins stay close to the mean.

files = pyalps.getResultFiles(prefix='parm1a')

# Check for linear drift in the measurement time series. Passing this check
# does not prove equilibration: autocorrelation and a stuck chain require
# separate inspection. A higher confidence level widens the acceptance region
# for zero slope; it is not a stricter equilibration check.
#
# Here is one example:
print(pyalps.checkSteadyState(outfile=files[0], observable='|Magnetization|', confidenceInterval=0.95))
print()

# and another one:
observables = pyalps.loadMeasurements(files, ['|Magnetization|', 'Energy'])
observables = pyalps.checkSteadyState(observables, confidenceInterval=0.95)
for o in observables:
    print('{}:\t{}'.format(o.props['observable'], o.props['checkSteadyState']))

# Check whether the error estimates have converged; if not, increase SWEEPS:
observables = pyalps.loadMeasurements(files, '|Magnetization|')
observables = pyalps.checkConvergence(observables)
for o in observables:
    print('{}:\t{}'.format(o.props['observable'], o.props['checkConvergence']))
