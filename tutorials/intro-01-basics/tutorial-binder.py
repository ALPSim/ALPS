from __future__ import print_function
# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Matthias Troyer <troyer@phys.ethz.ch> 
# 
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
import numpy as np

#prepare the input parameters
parms = []
for t in [1.5,2,2.5]:
    parms.append(
        { 
          'LATTICE'        : "square lattice", 
          'T'              : t,
          'J'              : 1 ,
          'THERMALIZATION' : 1000,
          'SWEEPS'         : 100000,
          'UPDATE'         : "cluster",
          'MODEL'          : "Ising",
          'L'              : 8
        }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm1',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5,writexml=True)

#get the list of result files
result_files = pyalps.getResultFiles(prefix='parm1')
print("Loading results from the files: ", result_files)

#print the observables stored in those files:
print("The files contain the following mesurements:", end=' ')
print(pyalps.loadObservableList(result_files))

#load a selection of measurements:
data = pyalps.loadMeasurements(result_files,['|Magnetization|','Magnetization^2'])

obschoose = lambda d, o: np.array(d)[np.nonzero([xx.props['observable'] == o for xx in d])]
binder = []
for dd in data:
    magn2 = obschoose(dd, 'Magnetization^2')[0]
    magnabs = obschoose(dd, '|Magnetization|')[0]
    
    res = pyalps.DataSet()
    res.props = pyalps.dict_intersect([d.props for d in dd])
    res.x = np.array([magnabs.props['T']])
    res.y = np.array([magn2.y[0]/(magnabs.y[0]*magnabs.y[0])])
    res.props['observable'] = 'Binder cumulant'
    binder.append(res)

binder = pyalps.collectXY(binder, 'T', 'Binder cumulant')

# ... and plot them
plt.figure()
pyalps.plot.plot(binder)
plt.xlabel('T')
plt.ylabel('Binder cumulant')
plt.show()
