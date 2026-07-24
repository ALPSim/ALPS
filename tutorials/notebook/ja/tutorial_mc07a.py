# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Brigitte Surer <surerb@phys.ethz.ch> 
# 
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import matplotlib.pyplot as plt
import pyalps.plot

#prepare the input parameters
parms = []
for l in [4,8,16]: 
    for t in [5.0,4.5,4.0,3.5,3.0,2.9,2.8,2.7]:
        parms.append(
            { 
              'LATTICE'        : "square lattice", 
              'T'              : t,
              'J'              : 1 ,
              'THERMALIZATION' : 1000,
              'SWEEPS'         : 400000,
              'UPDATE'         : "cluster",
              'MODEL'          : "Ising",
              'L'              : l
            }
    )
    for t in [2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.2]:
        parms.append(
            { 
              'LATTICE'        : "square lattice", 
              'T'              : t,
              'J'              : 1 ,
              'THERMALIZATION' : 1000,
              'SWEEPS'         : 40000,
              'UPDATE'         : "cluster",
              'MODEL'          : "Ising",
              'L'              : l
            }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('mc07a',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)
# use the following instead if you have MPI
#pyalps.runApplication('spinmc',input_file,Tmin=5,MPI=2)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='mc07a'))

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='mc07a'),['|Magnetization|', 'Connected Susceptibility', 'Specific Heat', 'Binder Cumulant', 'Binder Cumulant U2'])
magnetization_abs = pyalps.collectXY(data,x='T',y='|Magnetization|',foreach=['L'])
connected_susc = pyalps.collectXY(data,x='T',y='Connected Susceptibility',foreach=['L'])
spec_heat = pyalps.collectXY(data,x='T',y='Specific Heat',foreach=['L'])
binder_u4 = pyalps.collectXY(data,x='T',y='Binder Cumulant',foreach=['L'])
binder_u2 = pyalps.collectXY(data,x='T',y='Binder Cumulant U2',foreach=['L'])

#make plots
plt.figure()
pyalps.plot.plot(magnetization_abs)
plt.xlabel('Temperature $T$')
plt.ylabel('Magnetization $|m|$')
plt.title('2D Ising model')

plt.figure()
pyalps.plot.plot(connected_susc)
plt.xlabel('Temperature $T$')
plt.ylabel('Connected Susceptibility $\chi_c$')
plt.title('2D Ising model')

plt.figure()
pyalps.plot.plot(spec_heat)
plt.xlabel('Temperature $T$')
plt.ylabel('Specific Heat $c_v$')
plt.title('2D Ising model')

plt.figure()
pyalps.plot.plot(binder_u4)
plt.xlabel('Temperature $T$')
plt.ylabel('Binder Cumulant U4 $g$')
plt.title('2D Ising model')

plt.figure()
pyalps.plot.plot(binder_u2)
plt.xlabel('Temperature $T$')
plt.ylabel('Binder Cumulant U2 $g$')
plt.title('2D Ising model')
plt.show()

