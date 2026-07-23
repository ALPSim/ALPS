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

#get the list of result files
result_files = pyalps.getResultFiles(prefix='parm1')
print("Loading results from the files: ", result_files)

#print the observables stored in those files:
print("The files contain the following mesurements:", end=' ')
print(pyalps.loadObservableList(result_files))

#load a selection of measurements:
data = pyalps.loadMeasurements(result_files,['|Magnetization|','Magnetization^2'])

#make a plot for the magnetization: collect Magnetziation as function of T
plotdata = pyalps.collectXY(data,'T','|Magnetization|')
plt.figure()
pyalps.plot.plot(plotdata)
plt.xlim(0,3)
plt.ylim(0,1)
plt.title('Ising model')
plt.show()
