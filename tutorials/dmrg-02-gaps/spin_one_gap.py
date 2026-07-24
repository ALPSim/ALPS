# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2025 by ALPS Collaboration
#
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import numpy as np

#prepare the input parameters
parms = [ { 
        'LATTICE'                   : "open chain lattice", 
        'MODEL'                     : "spin",
        'local_S'                   : '1',
        'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
        'Sz_total'                  : 0,
        'J'                         : 1,
        'SWEEPS'                    : 5,
        'L'                         : 64,
        'MAXSTATES'                 : 300,
        'NUMBER_EIGENVALUES'        : 4
       } ]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm_spin_one_gap',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm_spin_one_gap'))

# print results:
energies = np.empty(0)
for s in data[0]:
    if s.props['observable'] == 'Energy':
        energies = s.y
    else:
        print(s.props['observable'], ':', s.y[0])
energies.sort()
print('Energies:', end=' ')
for e in energies:
    print(e, end=' ')
print('\nGap:', abs(energies[1]-energies[0]), abs(energies[2]-energies[1]))
