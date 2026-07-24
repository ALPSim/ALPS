from __future__ import print_function
#############################################################################
#
# ALPS Project Applications: Directed Worm Algorithm  
#
# Copyright (C) 2013 by Lode Pollet      <pollet@phys.ethz.ch>  
#                       Ping Nang Ma     <pingnang@phys.ethz.ch> 
#                       Matthias Troyer  <troyer@phys.ethz.ch>    
#
# SPDX-License-Identifier: MIT
#
#############################################################################

# The headers
import pyalps

# Set up a python list of parameters (python) dictionaries:
parms = [{
  'LATTICE'         : "square lattice",          
  'MODEL'           : "Ising",
  'L'               : 48,
  'J'               : 1.,
  'T'               : 2.269186,
  'THERMALIZATION'  : 10000,
  'SWEEPS'          : 50000,
}]

# Write into XML input file:
input_file = pyalps.writeInputFiles('parm1a',parms)

# and run the application spinmc:
pyalps.runApplication('spinmc', input_file, Tmin=10, writexml=True)

# We first get the list of all hdf5 result files via:
files = pyalps.getResultFiles(prefix='parm1a')

# and then extract, say the timeseries of the |Magnetization| measurements:
ts_M = pyalps.loadTimeSeries(files[0], '|Magnetization|');

# We can then visualize graphically:
import matplotlib.pyplot as plt
plt.plot(ts_M)
plt.show()

# ALPS Python provides a convenient tool to check whether a measurement observable(s) has (have) reached steady state equilibrium.
#
# Here is one example:
print(pyalps.checkSteadyState(outfile=files[0], observable='|Magnetization|', confidenceInterval=0.95))
print()

# and another one:
observables = pyalps.loadMeasurements(files, ['|Magnetization|', 'Energy'])
observables = pyalps.checkSteadyState(observables, confidenceInterval=0.95)
for o in observables:
    print('{}:\t{}'.format(o.props['observable'], o.props['checkSteadyState']))

