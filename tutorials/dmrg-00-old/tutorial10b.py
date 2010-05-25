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
parms = []
for h in [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]:
    parms.append(
        { 
          'LATTICE'                        : "open chain lattice", 
          'MODEL'                          : "spin",
          'CONSERVED_QUANTUMNUMBERS'       : 'Sz',
          'SWEEPS'                         : 4,
          'MAXSTATES'                      : 100,
          'L'                              : 10,
          'local_S'                        : 0.5,
          'J'                              : 1,
          'h'                              : h,
          'MEASURE_AVERAGE[Magnetization]' : 'Sz'
        }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm10b',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load the magnetization and collect it as function of field h
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm10b'))
magnetization = pyalps.collectXY(data,x='h',y='Magnetization')

magnetization[0].y /= magnetization[0].props['L']
#make plot
plt.figure()
pyalps.pyplot.plot(magnetization)
plt.xlabel('Field $h$')
plt.ylabel('Magnetization density $m$')
plt.ylim(0.0,0.5)
plt.title('Quantum Heisenberg chain')
