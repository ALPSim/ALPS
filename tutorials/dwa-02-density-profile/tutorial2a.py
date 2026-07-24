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

parms = [
  {
    'LATTICE' : 'inhomogeneous simple cubic lattice' ,
    'L'       : 21 ,

    'MODEL'   : 'boson Hubbard' ,
    'Nmax'    : 5 ,

    't'  : 1. ,
    'U'  : 8.11 ,
    'mu' : '4.05 - (0.2*(x-(L-1)/2.)*(x-(L-1)/2.) + 0.2*(y-(L-1)/2.)*(y-(L-1)/2.) + 0.2*(z-(L-1)/2.)*(z-(L-1)/2.))' ,

    'T'  : 1. ,

    'THERMALIZATION' : 2000,
    'SWEEPS'         : 20000,
    'SKIP'           : 10, 

    'MEASURE[Local Density]': 1
  }
]

input_file = pyalps.writeInputFiles('parm2a', parms)
res = pyalps.runApplication('dwa', input_file)


# Evaluating and plotting in Python
import pyalps
import pyalps.plot as aplt;

data     = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm2a'), 'Local Density');
aplt.plot3D(data, centeredAtOrigin=True)

