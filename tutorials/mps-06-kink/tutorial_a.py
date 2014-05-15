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

basename = 'sim_a'

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
for tau in [20.0]:
        ns = 500
        dt = tau / ns
        p = deepcopy(model)
        p['init_state'      ] = 'local_quantumnumbers'
        p['initial_local_S' ] = ','.join(['0.5']*50)
        p['initial_local_Sz'] = ','.join(['-0.5']*25 + ['0.5']*25)
        p['te_order' ] = 'second'
        p['dt'       ] = dt
        p['SWEEPS'   ] = ns
        p['tau'      ] = tau # not used in the simulation, but useful in the evaluation below
        p['always_measure'] = 'Local Magnetization'
        p['chkp_each'     ] = ns
        p['measure_each'  ] = 5
        p['COMPLEX'       ] = 1
        
        parms.append(p)


## write input files and run application
input_file = pyalps.writeInputFiles(basename, parms)
res = pyalps.runApplication('mps_evolve', input_file)


## simulation results
data = pyalps.loadIterationMeasurements(pyalps.getResultFiles(prefix=basename), what=['Local Magnetization'])

numeric_magnetization = []
for d in pyalps.flatten(data):
    L = d.props['L']
    for loc in [1,2]:
        q = deepcopy(d)
        q.x = [0]
        q.y = np.array([ q.y[0][L/2-loc] ])
        q.props['loc'] = loc
        numeric_magnetization.append(q)

mag_vs_time = pyalps.collectXY(numeric_magnetization, x='iter_t', y='Local Magnetization', foreach=['loc'])
for d in mag_vs_time:
    d.x = (d.x + 1.) * d.props['dt'] # convert time index to real time
    d.props['label'] = 'Numerical at n='+str(d.props['loc'])

plt.figure()
pyalps.plot.plot(mag_vs_time)
plt.xlabel('Time $t$')
plt.ylabel('Magnetization')
plt.title('Magnetization vs. time')
plt.legend(loc='best')

plt.show()
