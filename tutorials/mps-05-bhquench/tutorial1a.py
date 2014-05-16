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

basename = 'tutorial_1a'

## common model parameters
model = OrderedDict()
model['LATTICE'] = 'open chain lattice'
model['L'      ] = 10
model['MODEL'                   ] = 'hardcore boson'
model['CONSERVED_QUANTUMNUMBERS'] = 'N'
model['N_total'                 ] = 5
model['t'                       ] = 1.
model['MAXSTATES'] = 40

## ground state simulation for ground state
ground_parms = deepcopy(model)
ground_parms['V'] = 10.
ground_parms['SWEEPS'] = 4
input_file = pyalps.writeInputFiles(basename+'.ground',[ground_parms])
res = pyalps.runApplication('mps_optim',input_file)

initstate = pyalps.getResultFiles(prefix=basename+'.ground')[0].replace('xml', 'chkp')


def quench(g_i, g_f, steps, dt, p=1):
    t = ( 1. + np.arange(steps) ) * dt
    return g_i + (t/t[-1])**p * (g_f - g_i)
def list_to_string(ll):
    return ','.join([ str(li) for li in ll ])

#prepare the input parameters
parms = []
for A in [5.0, 10.0, 15.0, 25.0, 50.0]:
        ns = 500 + 500
        dt = 2*A / ns
        V = np.concatenate( [quench(10., 0.,  500, dt),   # first part quench from V=10 to V=0
                             quench(0.,  10., 500, dt)] ) # second part quench from V=0 to V=10
        
        p = deepcopy(model)
        p['initfile'] = initstate
        p['te_order' ] = 'second'
        p['DT'       ] = dt
        p['TIMESTEPS'] = ns
        p['V[Time]'  ] = list_to_string(V)
        p['tau'      ] = A # not used in the simulation, but useful in the evaluation below
        p['MEASURE_OVERLAP[Overlap]'] = initstate
        p['always_measure'] = 'Overlap'
        p['chkp_each'     ] = ns
        p['measure_each'  ] = 10
        p['update_each'   ] = 1
        p['COMPLEX'       ] = 1
        
        parms.append(p)


## write input files and run application
input_file = pyalps.writeInputFiles(basename+'.dynamic', parms)
res = pyalps.runApplication('mps_evolve', input_file)


## simulation results
data = pyalps.loadIterationMeasurements(pyalps.getResultFiles(prefix=basename+'.dynamic'), what=['Overlap'])


LE = pyalps.collectXY(data, x='Time', y='Overlap', foreach=['tau'])
for d in pyalps.flatten(LE):
    d.x =  (d.x + 1.) * d.props['dt'] # convert time index to real time
    d.y = abs(d.y)**2 # Loschmidt Echo defined as the module squared of the overlap
    d.props['label']=r'$\tau={0}$'.format( d.props['tau'] )

plt.figure()
pyalps.plot.plot(LE)
plt.xlabel('Time $t$')
plt.ylabel('Loschmidt Echo $|< \psi(0)|\psi(t) > |^2$')
plt.title('Loschmidt Echo vs. Time')
plt.legend(loc='lower right')


## Read V[Time] from props
Ufig = pyalps.collectXY(data, x='Time', y='V', foreach=['tau'])
for d in pyalps.flatten(Ufig):
    d.x =  (d.x + 1.) * d.props['dt'] # convert time index to real time
    d.props['label']=r'$\tau={0}$'.format( d.props['tau'] )

plt.figure()
pyalps.plot.plot(Ufig)
plt.xlabel('Time $t$')
plt.ylabel('V')
plt.title('Interaction parameter $V$ vs. Time')
plt.legend(loc='lower right')

plt.show()

