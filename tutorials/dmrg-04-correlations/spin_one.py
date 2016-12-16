from __future__ import division
# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Jan Gukelberger <gukelberger@phys.ethz.ch> 
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
parms = []
L = 32
for D in [20,40,60]:
    parms.append( { 
        'LATTICE_LIBRARY'                       : 'my_lattices.xml',
        'LATTICE'                               : 'open chain lattice with special edges '+str(L),
        'MODEL'                                 : 'spin',
        'local_S0'                              : 0.5,
        'local_S1'                              : 1,
        'CONSERVED_QUANTUMNUMBERS'              : 'N,Sz',
        'Sz_total'                              : 0,
        'J'                                     : 1,
        'SWEEPS'                                : 4,
        'NUMBER_EIGENVALUES'                    : 1,
        'MAXSTATES'                             : D,
        'MEASURE_AVERAGE[Magnetization]'        : 'Sz',
        'MEASURE_AVERAGE[Exchange]'             : 'exchange',
        'MEASURE_LOCAL[Local magnetization]'    : 'Sz',
        'MEASURE_CORRELATIONS[Diagonal spin correlations]'      : 'Sz',
        'MEASURE_CORRELATIONS[Offdiagonal spin correlations]'   : 'Splus:Sminus'
       } )


#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm_spin_one',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm_spin_one'))

# extract Sz correlation data
curves = []
for run in data:
    for s in run:
        if s.props['observable'] == 'Diagonal spin correlations':
            d = pyalps.DataSet()
            d.props['observable'] = 'Sz correlations'
            d.props['label'] = 'D = '+str(s.props['MAXSTATES'])
            d.x = np.arange(L)
            
            # sites with increasing distance l symmetric to the chain center
            site1 = np.array([int(-(l+1)/2.0) for l in range(0,L)]) + L//2
            site2 = np.array([int(  l   /2.0) for l in range(0,L)]) + L//2
            indices = L*site1 + site2
            d.y = abs(s.y[0][indices])
            
            curves.append(d)

# Plot correlation vs. distance
plt.figure()
pyalps.plot.plot(curves)
plt.yscale('log')
plt.legend()
plt.title('Spin correlations in antiferromagnetic Heisenberg chain (S=1)')
plt.ylabel('correlations $| \\langle S^z_{L/2-l/2} S^z_{L/2+l/2} \\rangle |$')
plt.xlabel('distance $l$')

plt.show()
