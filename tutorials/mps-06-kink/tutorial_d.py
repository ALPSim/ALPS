# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2014 by Michele Dolfi <dolfim@phys.ethz.ch>
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
import pyalps.plot
import numpy as np
from collections import OrderedDict
from copy import deepcopy

basename = 'sim_d'

## common model parameters
model = OrderedDict()
model['LATTICE'] = 'open chain lattice'
model['L'      ] = 50
model['MODEL'                   ] = 'spin'
model['CONSERVED_QUANTUMNUMBERS'] = 'Sz'
model['Sz_total'                ] = 0
model['Jxy'                     ] = 1.
model['MAXSTATES'] = 40
model['MEASURE_LOCAL[Local Magnetization]'] = 'Sz'

#prepare the input parameters
parms = []
for z in [0.0, 0.3, 0.9, 1.0, 1.1, 1.5]:
    tau = 20.
    nsteps = 500
    dt = tau / nsteps
    p = deepcopy(model)
    p['init_state'      ] = 'local_quantumnumbers'
    p['initial_local_S' ] = ','.join(['0.5']*50)
    p['initial_local_Sz'] = ','.join(['-0.5']*25 + ['0.5']*25)
    p['te_order' ] = 'second'
    p['DT'       ] = dt
    p['TIMESTEPS'] = nsteps
    p['tau'      ] = tau # not used in the simulation, but useful in the evaluation below
    p['Jz'       ] = z
    p['always_measure'] = 'Local Magnetization'
    p['chkp_each'     ] = nsteps
    p['measure_each'  ] = 5
    p['COMPLEX'       ] = 1
    
    parms.append(p)


## write input files and run application
input_file = pyalps.writeInputFiles(basename, parms)
res = pyalps.runApplication('mps_evolve', input_file)


## simulation results
data = pyalps.loadIterationMeasurements(pyalps.getResultFiles(prefix=basename), what=['Local Magnetization'])

for q in pyalps.flatten(data):
    L=q.props['L']
    #Compute the integrated flow of magnetization through the center \Delta M=\sum_{n>L/2}^{L} (<S_n^z(t)>+1/2)
    #\Delta M= L/4
    loc=0.5*(L/2)
    #\Delta M-=<S_n^z(t)> from n=L/2 to L
    q.y=np.array([0.5*(L/2)-sum(q.y[0][L/2:L])])


#Plot the Error in the magnetization one site to the right of the chain center
Mag=pyalps.collectXY(data, x='Time', y='Local Magnetization', foreach=['Jz'])
for d in Mag:
    d.x = (d.x+1) * d.props['DT']

plt.figure()
pyalps.plot.plot(Mag)
plt.xlabel('Time $t$')
plt.ylabel('Integrated Magnetization $\Delta M(t)$')
plt.title('Integrated Magnetization vs. Time')
plt.legend(loc='upper left')
plt.show()

