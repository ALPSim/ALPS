# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Jan Gukelberger <gukelberger@phys.ethz.ch> 
# 
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.plot

#prepare the input parameters
parms = []
for sz in [0,1,2]:
    parms.append( { 
        'LATTICE_LIBRARY'           : 'my_lattices.xml',
        'LATTICE'                   : 'open chain lattice with special edges 32',
        'MODEL'                     : "spin",
        'local_S0'                  : '0.5',
        'local_S1'                  : '1',
        'CONSERVED_QUANTUMNUMBERS'  : 'N,Sz',
        'Sz_total'                  : sz,
        'J'                         : 1,
        'SWEEPS'                    : 4,
        'NUMBER_EIGENVALUES'        : 1,
        'MAXSTATES'                 : 40,
        'MEASURE_LOCAL[Local magnetization]'   : 'Sz'
    } )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm_spin_one',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm_spin_one'))

# extract local magnetization data
curves = []
for run in data:
    for s in run:
        if s.props['observable'] == 'Local magnetization':
            sz = s.props['Sz_total']
            s.props['label'] = '$S_z = ' + str(sz) + '$'
            s.y = s.y.flatten()
            curves.append(s)
        
# Plot local magnetization vs. site
plt.figure()
pyalps.plot.plot(curves)
plt.legend()
plt.title('Magnetization of antiferromagnetic Heisenberg chain (S=1)')
plt.ylabel('local magnetization')
plt.xlabel('site')

plt.show()
