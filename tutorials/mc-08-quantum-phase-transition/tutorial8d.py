# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Brigitte Surer
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
import numpy as np

#prepare the input parameters
parms = []
for l in [16,32,64,128]:
    # the two largest sizes only scan the immediate neighbourhood of the crossing
    if l <= 32:
        j2list = [0.31,0.311,0.312,0.313,0.314,0.315,0.316,0.317,0.318,0.319,0.32]
    else:
        j2list = [0.312,0.313,0.314,0.315,0.316]
    for j2 in j2list:
        parms.append(
            {
              'LATTICE'        : "coupled ladders",
              'local_S'        : 0.5,
              'ALGORITHM'      : 'loop',
              'SEED'           : 0,
              'BETA'           : 2*l,
              'J0'             : 1 ,
              'J1'             : 1,
              'J2'             : j2,
              'THERMALIZATION' : 5000,
              'SWEEPS'         : 50000,
              'MODEL'          : "spin",
              'L'              : l,
              'W'              : l//2
            }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm8d',parms)
pyalps.runApplication('loop',input_file)

data = pyalps.loadMeasurements(pyalps.getResultFiles(pattern='parm8d.task*.out.h5'),['Binder Ratio of Staggered Magnetization','Stiffness'])

binder=pyalps.collectXY(data,x='J2',y='Binder Ratio of Staggered Magnetization', foreach=['L'])
stiffness =pyalps.collectXY(data,x='J2',y='Stiffness', foreach=['L'])

for q in stiffness:
    q.y = q.y*q.props['L']

#make plot    
plt.figure()
pyalps.plot.plot(stiffness)
plt.xlabel(r'$J2$')
plt.ylabel(r'Stiffness $\rho_s L$')
plt.title('coupled ladders')

plt.figure()
pyalps.plot.plot(binder)
plt.xlabel(r'$J_2$')
plt.ylabel(r'$g(m_s)$')
plt.title('coupled ladders')
plt.show()
