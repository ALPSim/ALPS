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
for D in range(50,250,50):
    for L in range(32,160,32):
        for sz in [0,1]:
            parms.append( { 
        'LATTICE'                   : "open chain lattice", 
        'MODEL'                     : "spin",
        'CONSERVED_QUANTUMNUMBERS'  : 'N,Sz',
        'Sz_total'                  : sz,
        'J'                         : 1,
        'SWEEPS'                    : 4,
        'NUMBER_EIGENVALUES'        : 1,
        'L'                         : L,
        'MAXSTATES'                 : D,
        'NUMBER_EIGENVALUES'        : 1
       } )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm_spin_one_half_triplet',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm_spin_one_half_triplet'))

# extract energies
dvals = []
lvals = []
energies = {}
for run in data:
    for s in run:
        if s.props['observable'] == 'Energy':
            D = s.props['MAXSTATES']
            L = s.props['L']
            sz = s.props['Sz_total']
            dvals.append(D)
            lvals.append(L)
            energies[(D,L,sz)] = s.y[0]
        
# Plot gap vs. 1/L for fixed D
lplot = []
for D in dvals:
    curve = pyalps.DataSet()
    curve.props['label'] = 'D = ' + str(D)
    curve.x = np.array([1/L for L in lvals])
    curve.y = np.array([energies[(D,L,1)]-energies[(D,L,0)] for L in lvals])
    lplot.append(curve)
plt.figure()
pyalps.pyplot.plot(lplot)
plt.legend()
plt.title('Gap of antiferromagnetic Heisenberg chain (S=1/2)')
plt.ylabel('Gap $\Delta$')
plt.xlabel('$1/L$')
#plt.xlim(0,2*3.1416)
#plt.ylim(0,2.5)



plt.show()
