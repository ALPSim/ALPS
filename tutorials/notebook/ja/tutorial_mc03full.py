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

# Please run the two other tutorials before running this one. 
# This tutorial relies on the results created in those tutorials

import pyalps
import matplotlib.pyplot as plt
import pyalps.plot

# load all files
data = pyalps.loadMeasurements(pyalps.getResultFiles(),'Magnetization Density')

#flatten the hierarchical structure
data = pyalps.flatten(data)

#load the magnetization and collect it as function of field h
magnetization = pyalps.collectXY(data,x='h',y='Magnetization Density',foreach=['LATTICE'])

#make plot
plt.figure()
pyalps.plot.plot(magnetization)
plt.xlabel('Field $h$')
plt.ylabel('Magnetization $m$')
plt.ylim(0.0,0.5)
plt.legend()
plt.show()
