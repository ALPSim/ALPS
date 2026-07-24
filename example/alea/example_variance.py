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

# This is an example of how to easily calculate the variance of data stored in a hdf5 file.


filename = "testfile.h5"

# create the correct MCData object to load the data.
obs = pyalps.alea.MCScalarData()

# load the variable E saved in the file testfile.h5 into the mcdata object.
obs.load(filename, "simulation/results/" + pyalps.hdf5_name_encode("E"))

# calculate the mean
variance = pyalps.alea.variance(obs)

# print the result
print("The variance of E is: " + str(variance))
  
# write the result back to the file
ar = pyalps.hdf5.archive(filename, 1)
ar.write("simulation/results/" + pyalps.hdf5_name_encode("E") + "/variance/value", variance)


