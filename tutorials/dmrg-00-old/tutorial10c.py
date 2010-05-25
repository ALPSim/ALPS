# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Matthias Troyer <troyer@phys.ethz.ch> 
# 
# This software is part of the ALPS libraries, published under the ALPS
# Library License; you can use, redistribute it and/or modify it under
# the terms of the license, either version 1 or (at your option) any later
# version.
#  
# You should have received a copy of the ALPS Library License along with
# the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
# available from http://alps.comp-phys.org/.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
# FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
# 
# ****************************************************************************

import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot
import numpy as np

#prepare the input parameters
parms = [{ 
          'LATTICE'                   : "open chain lattice", 
          'MODEL'                     : "spin",
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'SWEEPS'                    : 4,
          'MAXSTATES'                 : 100,
          'L'                         : 10,
          'local_S'                   : 1,
          'Sz_total'                  : 0,
          'J'                         : 1,
          'h'                         : 0,
          'NUMBER_EIGENVALUES'        : 2,
          'MEASURE_AVERAGE[Magnetization]'                      : 'Sz',
          'MEASURE_AVERAGE[Exchange]'                           : 'exchange',
          'MEASURE_LOCAL[Local magnetization]'                  : 'Sz',
          'MEASURE_CORRELATIONS[Diagonal spin correlations]='   : 'Sz',
          'MEASURE_CORRELATIONS[Offdiagonal spin correlations]' : 'Splus:Sminus'

        }]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm10c',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm10c'))

# print properties of the ground state:
for s in data[0]:
  if pyalps.size(s.y[0])==1:
    print s.props['observable'], ' : ', s.y[0]
  else:
    for (x,y) in zip(s.x,s.y[0]):
      print  s.props['observable'], '(', x, ') : ', y
