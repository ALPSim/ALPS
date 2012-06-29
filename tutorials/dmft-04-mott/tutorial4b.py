# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Brigitte Surer <surerb@phys.ethz.ch> 
#               2012 by Jakub Imriska  <jimriska@phys.ethz.ch>
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
for u in [4.,5.,6.,8.]: 
  for b in [20.]:
    parms.append(
            { 
              'ANTIFERROMAGNET'         : 0,
              'CHECKPOINT'              : 'solverdump_U_'+str(u),
              'CONVERGED'               : 0.001,
              'FLAVORS'                 : 2,
              'H'                       : 0,
              'H_INIT'                  : 0.,
              'MAX_IT'                  : 20,
              'MAX_TIME'                : 10,
              'MU'                      : 0,
              'N'                       : 500,
              'NMATSUBARA'              : 500, 
              'N_MEAS'                  : 10000,
              'N_ORDER'                 : 50,
              'OMEGA_LOOP'              : 1,
              'SEED'                    : 0, 
              'SITES'                   : 1,              
              'SOLVER'                  : 'Hybridization',
              'SYMMETRIZATION'          : 1,
              't'                       : 1,
              'SWEEPS'                  : 100000000,
              'BETA'                    : b,
              'THERMALIZATION'          : 1000,
              'U'                       : u
            }
        )

# NOTE: in revision of ALPS older than 6238, the MAX_TIME will effectively be 60 seconds.        
# For more precise calculations change the parameters
#   MAX_TIME (to 60-180),
#   CONVERGED (to 0.0002-0.001),
#   raise N and NMATSUBARA (to 1000)

## Please run the tutorial4a.py before this one or uncomment the following lines.
## This tutorial relies on the results created there.

# #write the input file and run the simulation
# for p in parms:
#     input_file = pyalps.writeParameterFile('parm_u_'+str(p['U']),p)
#     res = pyalps.runDMFT(input_file)

listobs=['Green_0']   # we look at convergence of a single flavor (=0)
    
ll=pyalps.load.Hdf5Loader()
for p in parms:
    data = ll.ReadDMFTIterations(pyalps.getResultFiles(pattern='parm_u_'+str(p['U'])+'.h5'), measurements=listobs, verbose=True)
    grouped = pyalps.groupSets(pyalps.flatten(data), ['iteration'])
    nd=[]
    for group in grouped:
        r = pyalps.DataSet()
        r.y = -np.array(group[0].y)
        r.x = np.array([e*group[0].props['BETA']/float(group[0].props['N']) for e in group[0].x])
        r.props = group[0].props
        r.props['label'] = 'it'+r.props['iteration']
        nd.append( r )
    plt.figure()
    plt.yscale('log')
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$G_{flavor=0}(\tau)$')
    plt.title('DMFT-04: ' + r'$U = %.4s$' %nd[0].props['U'])
    pyalps.plot.plot(nd)
    plt.legend()
plt.show()

