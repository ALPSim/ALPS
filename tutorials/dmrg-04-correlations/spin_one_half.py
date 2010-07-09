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
import pyalps.pyplot

#prepare the input parameters
parms = []
for D in [10,30,50]: #[100,150]: #,200,250,300,400]:
    parms.append( { 
        'LATTICE'                               : 'open chain lattice', 
        'MODEL'                                 : 'spin',
        'CONSERVED_QUANTUMNUMBERS'              : 'N,Sz',
        'Sz_total'                              : 0,
        'J'                                     : 1,
        'SWEEPS'                                : 4, #6
        'NUMBER_EIGENVALUES'                    : 1,
        'L'                                     : 64, #192,
        'MAXSTATES'                             : D,
        'MEASURE_AVERAGE[Magnetization]'        : 'Sz',
        'MEASURE_AVERAGE[Exchange]'             : 'exchange',
        'MEASURE_LOCAL[Local magnetization]'    : 'Sz',
        'MEASURE_CORRELATIONS[Diagonal spin correlations]'      : 'Sz',
        'MEASURE_CORRELATIONS[Offdiagonal spin correlations]'   : 'Splus:Sminus'
       } )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm_spin_one_half',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm_spin_one_half'))

# extract Sz correlation data
curves = []
for run in data:
    for s in run:
        if s.props['observable'] == 'Diagonal spin correlations':
            d = pyalps.DataSet()
            d.props['observable'] = 'Sz correlations'
            d.props['label'] = 'D = '+str(s.props['MAXSTATES'])
            L = int(s.props['L'])
            print len(s.x), len(s.y)
            print s.y
            d.x = np.arange(L)
            d.y = np.array([s.y[0][L*(L/2+int(-(l+1)/2.0))+L/2+int(l/2.0)] for l in range(0,L)])
            d.y = abs(d.y)
            curves.append(d)
            
#            for l in range(0,L):
#                i = L/2+int(-(l+1)/2.0)
#                j = L/2+int(l/2.0)
#                x = i*L+j
#                print l, j-i, i, j, s.x[x]
#            sz = s.props['Sz_total']
#            s.props['label'] = '$S_z = ' + str(sz) + '$'
#            s.y = s.y.flatten()
#            curves.append(s)
        
# Plot correlation vs. distance
plt.figure()
pyalps.pyplot.plot(curves)
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.title('Spin correlations in antiferromagnetic Heisenberg chain (S=1/2)')
plt.ylabel('correlations $| \\langle S^z_{L/2-l/2} S^z_{L/2+l/2} \\rangle |$')
plt.xlabel('distance $l$')


plt.show()
