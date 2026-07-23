from __future__ import division
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
L = 32
for D in [20,40,60]:
    parms.append( { 
        'LATTICE_LIBRARY'                       : 'my_lattices.xml',
        'LATTICE'                               : 'open chain lattice with special edges '+str(L),
        'MODEL'                                 : 'spin',
        'local_S0'                              : 0.5,
        'local_S1'                              : 1,
        'CONSERVED_QUANTUMNUMBERS'              : 'N,Sz',
        'Sz_total'                              : 0,
        'J'                                     : 1,
        'SWEEPS'                                : 4,
        'NUMBER_EIGENVALUES'                    : 1,
        'MAXSTATES'                             : D,
        'MEASURE_AVERAGE[Magnetization]'        : 'Sz',
        'MEASURE_AVERAGE[Exchange]'             : 'exchange',
        'MEASURE_LOCAL[Local magnetization]'    : 'Sz',
        'MEASURE_CORRELATIONS[Diagonal spin correlations]'      : 'Sz',
        'MEASURE_CORRELATIONS[Offdiagonal spin correlations]'   : 'Splus:Sminus'
       } )


#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm_spin_one',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm_spin_one'))

# extract Sz correlation data
curves = []
for run in data:
    for s in run:
        if s.props['observable'] == 'Diagonal spin correlations':
            d = pyalps.DataSet()
            d.props['observable'] = 'Sz correlations'
            d.props['label'] = 'D = '+str(s.props['MAXSTATES'])
            d.x = np.arange(L)
            
            # sites with increasing distance l symmetric to the chain center
            site1 = np.array([int(-(l+1)/2.0) for l in range(0,L)]) + L//2
            site2 = np.array([int(  l   /2.0) for l in range(0,L)]) + L//2
            indices = L*site1 + site2
            d.y = abs(s.y[0][indices])
            
            curves.append(d)

# Plot correlation vs. distance
plt.figure()
pyalps.plot.plot(curves)
plt.yscale('log')
plt.legend()
plt.title('Spin correlations in antiferromagnetic Heisenberg chain (S=1)')
plt.ylabel('correlations $| \\langle S^z_{L/2-l/2} S^z_{L/2+l/2} \\rangle |$')
plt.xlabel('distance $l$')

plt.show()
