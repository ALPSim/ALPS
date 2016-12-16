from __future__ import print_function
# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Jan Gukelberger <gukelberger@phys.ethz.ch>
#               2014 by Michele Dolfi <dolfim@phys.ethz.ch>
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
import numpy as np
import matplotlib.pyplot as plt
import pyalps.plot

#prepare the input parameters
parms = [ { 
        'LATTICE'                   : "open chain lattice", 
        'MODEL'                     : "spin",
        'CONSERVED_QUANTUMNUMBERS'  : 'N,Sz',
        'Sz_total'                  : 0,
        'J'                         : 1,
        'SWEEPS'                    : 4,
        'NUMBER_EIGENVALUES'        : 1,
        'L'                         : 32,
        'MAXSTATES'                 : 100
       } ]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm_spin_one_half',parms)
res = pyalps.runApplication('mps_optim',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm_spin_one_half'))

# print properties of the eigenvector:
for s in data[0]:
    print(s.props['observable'], ' : ', s.y[0])

# load and plot iteration history
iterations = pyalps.loadIterationMeasurements(pyalps.getResultFiles(prefix='parm_spin_one_half'),
                                      what=['Energy','TruncatedWeight'])
energy_iteration     = pyalps.collectXY(pyalps.flatten(iterations), 'iteration', 'Energy')
for d in energy_iteration:
    d.x = range(0, len(d.y))
truncation_iteration = pyalps.collectXY(pyalps.flatten(iterations), 'iteration', 'TruncatedWeight')
for d in truncation_iteration:
    d.x = range(0, len(d.y))

plt.figure()
pyalps.plot.plot(energy_iteration)
plt.title('Iteration history of ground state energy (S=1)')
plt.ylabel('$E_0$')
plt.xlabel('iteration')

plt.figure()
pyalps.plot.plot(truncation_iteration)
plt.title('Iteration history of truncation error (S=1)')
plt.yscale('log')
plt.ylabel('error')
plt.xlabel('iteration')

plt.show()
