from __future__ import print_function
#/*****************************************************************************
#*
#* ALPS Project: Algorithms and Libraries for Physics Simulations
#*
#* Copyright (C) 2011-2012 by Lukas Gamper <gamperl@gmail.com>,
#*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
#*                            Maximilian Poprawe <poprawem@ethz.ch>
#*
#* SPDX-License-Identifier: MIT
#*
#*****************************************************************************/

import pyalps

# This is an example of how to easily calculate the mean of data stored in a hdf5 file.


filename = "testfile.h5"

# create the correct MCData object to load the data.
obs = pyalps.alea.MCScalarData()

# load the variable E saved in the file testfile.h5 into the mcdata object.
obs.load(filename, "simulation/results/" + pyalps.hdf5_name_encode("E"))

# calculate the mean
mean = pyalps.alea.mean(obs)

# print the result
print("The mean of E is: " + str(mean))
  
# write the result back to the file
ar = pyalps.hdf5.archive(filename, 1)
ar.write("simulation/results/" + pyalps.hdf5_name_encode("E") + "/mean/value", mean)


