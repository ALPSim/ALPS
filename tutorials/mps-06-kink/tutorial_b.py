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
import scipy.special

basename = 'sim_b'

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
for nsteps in [100, 250, 500, 750, 1000]:
    tau = 20.
    dt = tau / nsteps
    p = deepcopy(model)
    p['init_state'      ] = 'local_quantumnumbers'
    p['initial_local_S' ] = ','.join(['0.5']*50)
    p['initial_local_Sz'] = ','.join(['-0.5']*25 + ['0.5']*25)
    p['te_order' ] = 'second'
    p['DT'       ] = dt
    p['TIMESTEPS'] = nsteps
    p['tau'      ] = tau # not used in the simulation, but useful in the evaluation below
    p['always_measure'] = 'Local Magnetization'
    p['chkp_each'     ] = nsteps
    p['measure_each'  ] = 2
    p['MAXSTATES'     ] = 20
    p['COMPLEX'       ] = 1
    
    parms.append(p)


## write input files and run application
input_file = pyalps.writeInputFiles(basename, parms)
res = pyalps.runApplication('mps_evolve', input_file)


## simulation results
data = pyalps.loadIterationMeasurements(pyalps.getResultFiles(prefix=basename), what=['Local Magnetization'])

for q in pyalps.flatten(data):
    L=q.props['L']
    time=(q.props['Time']+1) * q.props['DT']
    #Get the exact result of M(1,t)=-(1/2)*(j_0(t)^2), where j_0(t) is the 0^{th} order
    # bessel function and M(1,t) is the magnetization one site to the right of the chain center
    loc=-0.5*scipy.special.jn(0,time)**2
    #Get the difference between the computed and exact results
    q.y=np.array([abs(q.y[0][L/2+1-1]-loc)])


#Plot the Error in the magnetization one site to the right of the chain center
Mag=pyalps.collectXY(data, x='Time', y='Local Magnetization', foreach=['DT'])
for d in Mag:
    d.x = (d.x+1) * d.props['DT']

plt.figure()
pyalps.plot.plot(Mag)
plt.xlabel('Time $t$')
plt.yscale('log')
plt.ylabel('Magnetization Error')
plt.title('Error in the magnetization vs. time')
plt.legend(loc='lower left')
plt.show()
