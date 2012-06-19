# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Brigitte Surer <surerb@phys.ethz.ch> 
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
parms=[]
coulombparam=[[1.8,0.45],[2.2,0.55],[2.8,0.7]]
for cp in coulombparam: 
    parms.append(
            { 
              'CHECKPOINT'          : 'dump',
              'CONVERGED'           : 0.01,
              'FLAVORS'             : 4,
              'H'                   : 0,
              'H_INIT'              : 0.,
              'MAX_IT'              : 20,
              'MAX_TIME'            : 180,
              'MU'                  : 0,
              'N'                   : 1000,
              'NMATSUBARA'          : 1000,
              'N_MEAS'              : 10000,
              'N_ORDER'             : 50,
              'OMEGA_LOOP'          : 1,
              'SEED'                : 0,
              'SOLVER'              : 'Hybridization',
              'SYMMETRIZATION'      : 1,
              'TOLERANCE'           : 0.3,
              't'                   : 1,
              'SWEEPS'              : 100000000,
              'BETA'                : 30,
              'THERMALIZATION'      : 10,
              'U'                   : cp[0],
              'J'                   : cp[1],
              't0'                  : 0.5,
              't1'                  : 1,
              'G0TAU_INPUT'         :'G0_tau_input_u_'+str(cp[0])+'_j_'+str(cp[1])
            }
        )

## Please run the tutorial5a.py before this one or uncomment the following lines.
## This tutorial relies on the results created there.

# #write the input file and run the simulation
# for p in parms:
#     input_file = pyalps.writeParameterFile('parm_u_'+str(p['U'])+'_j_'+str(p['J']),p)
#     res = pyalps.runDMFT(input_file)

flavors=parms[0]['FLAVORS']
listobs=[]   
for f in range(0,flavors):
    listobs.append('Green_'+str(f))
    
ll=pyalps.load.Hdf5Loader()
for p in parms:
    data = ll.ReadDMFTIterations(pyalps.getResultFiles(pattern='parm_u_'+str(p['U'])+'_j_'+str(p['J'])+'.h5'), measurements=listobs, verbose=True)
    grouped = pyalps.groupSets(pyalps.flatten(data), ['iteration'])
    nd=[]
    for group in grouped:
        r = pyalps.DataSet()
        r.y = np.array(group[0].y)
        r.x = np.array([e*group[0].props['BETA']/float(group[0].props['N']) for e in group[0].x])
        r.props = group[0].props
        r.props['label'] = r.props['iteration']
        nd.append( r )
    plt.figure()
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$G(\tau)$')
    plt.title(r'$U$ = %.4s' %nd[0].props['U'])
    pyalps.plot.plot(nd)
    plt.legend()

plt.show()
