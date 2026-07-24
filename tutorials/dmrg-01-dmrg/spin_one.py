from __future__ import print_function
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
parms = [ { 
        'LATTICE_LIBRARY'           : 'my_lattice.xml',
        'LATTICE'                   : 'open chain lattice with special edges',
        'MODEL'                     : 'spin',
        'local_S0'                  : '0.5',
        'local_S1'                  : '1',
        'CONSERVED_QUANTUMNUMBERS'  : 'N,Sz',
        'Sz_total'                  : 0,
        'J'                         : 1,
        'SWEEPS'                    : 4,
        'NUMBER_EIGENVALUES'        : 1,
        'MAXSTATES'                 : 100
       } ]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm_spin_one',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm_spin_one'))

# print properties of the eigenvector:
for s in data[0]:
    print(s.props['observable'], ' : ', s.y[0])

# load and plot iteration history
iter = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm_spin_one'),
                               what=['Iteration Energy','Iteration Truncation Error'])

plt.figure()
pyalps.plot.plot(iter[0][0])
plt.title('Iteration history of ground state energy (S=1)')
plt.ylabel('$E_0$')
plt.xlabel('iteration')

plt.figure()
pyalps.plot.plot(iter[0][1])
plt.title('Iteration history of truncation error (S=1)')
plt.yscale('log')
plt.ylabel('error')
plt.xlabel('iteration')

plt.show()
