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
import pyalps.plot

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

#make a plot for the magnetization: collect Magnetziation as function of T
plotdata = pyalps.collectXY(data,'T','|Magnetization|')

# convert the data to text file for plotting using another tool
print("The results in plain text are:")
print(pyalps.plot.convertToText(plotdata))

print("Saving into file parm1.txt")
f = open ('parm1.txt','w')
f.write(pyalps.plot.convertToText(plotdata))
f.close()
