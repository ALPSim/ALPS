from __future__ import print_function
#/*****************************************************************************
#*
#* ALPS Project: Algorithms and Libraries for Physics Simulations
#*
#* Copyright (C) 2011-2012 by Lukas Gamper <gamperl@gmail.com>,
#*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
#*                            Maximilian Poprawe <poprawem@ethz.ch>
#*
#* This software is part of the ALPS libraries, published under the ALPS
#* Library License; you can use, redistribute it and/or modify it under
#* the terms of the license, either version 1 or (at your option) any later
#* version.
#*
#* You should have received a copy of the ALPS Library License along with
#* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
#* available from http://alps.comp-phys.org/.
#*
#* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
#* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
#* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
#* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#* DEALINGS IN THE SOFTWARE.
#*
#*****************************************************************************/

import pyalps
import numpy

# This is an example of how to estimate the error of the mean of data stored in a hdf5 file

filename = "testfile.h5"

# create the correct MCData object to load the data.
obs = pyalps.alea.MCScalarData()

# load the variable m saved in the file testfile.h5 into the mcdata object.
obs.load(filename, "simulation/results/" + pyalps.hdf5_name_encode("m"))

# calculate the autocorrelation until it has reached 20% of its start value
auto_corr = pyalps.alea.autocorrelation(obs, _limit = 0.2)

# fit the autocorrelation exponentially between the values where it is at 80% and at 20% of the value at t = 1
fit = pyalps.alea.exponential_autocorrelation_time(auto_corr, _min=.2, _max=.8)

#calculate the integrated autocorrelation time by summing up the autocorrelation and integrating the fit for the tail
int_autocorr_time = pyalps.alea.integrated_autocorrelation_time(pyalps.alea.cut_tail(auto_corr, _limit=0.2), fit);

  # calculate different error estimates:
# assuming uncorrelated data:
error_uncorr = pyalps.alea.error(obs, pyalps.alea.uncorrelated)

# accounting for autocorrelation using binning analysis:
error_binning = pyalps.alea.error(obs, pyalps.alea.binning)

# accounting for autocorrelation using the integrated autocorrelation time:
error_corrtime = float(error_uncorr * numpy.sqrt(1. + 2. * int_autocorr_time))

# print the result
print("The estimated integrated autocorrelation time is: " + str(int_autocorr_time))
print("The different error estimates are:")
print("uncorrelated: " + str(error_uncorr))
print("with binning: " + str(error_binning))
print("with correlation time: " + str(error_corrtime))

# we can also write one of the errors back to the file
ar = pyalps.hdf5.archive(filename, 1)
ar.write("simulation/results/" + pyalps.hdf5_name_encode("m") + "/mean/error", error_corrtime)



