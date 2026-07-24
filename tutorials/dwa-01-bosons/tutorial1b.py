# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2013 by Matthias Troyer <troyer@phys.ethz.ch>,
#                       Ping Nang Ma    <pingnang@phys.ethz.ch> 
# 
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************


# Preparing and running the simulation using Python
import pyalps

parms = []
for L in [4,6,8]:
  for t in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
    parms.append(
        { 
          'LATTICE'                 : "square lattice", 
          'MODEL'                   : "boson Hubbard",
          'T'                       : 0.1,
          'L'                       : L ,
          't'                       : t ,
          'mu'                      : 0.5,
          'U'                       : 1.0 ,
          'Nmax'                    : 2 ,
          'THERMALIZATION'          : 100000,
          'SWEEPS'                  : 2000000,
          'SKIP'                    : 500,
          'MEASURE[Winding Number]': 1
        }
    )

input_file = pyalps.writeInputFiles('parm1b', parms)
res = pyalps.runApplication('dwa', input_file, Tmin=5, writexml=True)


# Evaluating the simulation and preparing plots using Python
import pyalps
import matplotlib.pyplot as plt
import pyalps.plot as aplt

data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm1b'),'Stiffness')
rhos = pyalps.collectXY(data,x='t',y='Stiffness',foreach=['L'])

for rho in rhos:
  rho.y = rho.y * float(rho.props['L'])

plt.figure()
aplt.plot(rhos)
plt.xlabel('Hopping $t/U$')
plt.ylabel('$\\rho _sL$')
plt.legend()
plt.title('Scaling plot for Bose-Hubbard model')
plt.show()

