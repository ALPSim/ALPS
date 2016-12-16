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


from optparse import OptionParser
from sys import exit
import pyalps
import pyalps.hdf5 as h5
import pyalps.alea as alea


def impl_calculation(name, save_path, calculate):

  usage = "Usage: %prog [options] FILE [FILE [...]]"
  parser = OptionParser(usage=usage)

  parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="print detailed information")
  parser.add_option("-w", "--write", action="store_true", dest="write", help="write the result(s) back into the file(s)")
  parser.add_option("-n", "--name", action="append", metavar="VAR", dest="variables", help="variable name, can be specified multiple times [default: all variables]")
  parser.add_option("-p", "--path", action="store", metavar="HDF5-PATH", dest="path", help="hdf5-path where the data is stored [default: \"/simulation/results\"]")
  parser.set_defaults(verbose = False, write = False, variables = [], path="/simulation/results")

  (options, args) = parser.parse_args()

  if len(args) == 0:
    parser.print_help()
    exit()

  variables = options.variables

  for filestring in args:
    ar = h5.archive(filestring, 1)
    if len(options.variables) == 0:
      variables = ar.list_children(options.path)
      if options.verbose:
        print("Variables in file " + filestring + ":  " + " ".join(variables))

    for variablestring in variables:
      if ar.dimensions(options.path + "/" + pyalps.hdf5_name_encode(variablestring) + "/timeseries/data") == 1:
        obs = alea.MCScalarData()
        #_save = mcanalyze.write_dim_0
        #vector_save = mcanalyze.write_dim_1
      else:
        obs = alea.MCVectorData()
        #scalar_save = mcanalyze.write_dim_1
        #vector_save = mcanalyze.write_dim_2
      obs.load(filestring, options.path +"/" + pyalps.hdf5_name_encode(variablestring))
      result = calculate(obs)
      if options.verbose:
        print("The " + name + " of variable " + variablestring + " in file " + filestring + " is: " + str(result))
      if options.write:
        ar.write(options.path + "/" + pyalps.hdf5_name_encode(variablestring) + "/" + save_path, result)
  print("Done")



