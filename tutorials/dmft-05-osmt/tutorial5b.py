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
coulombparam=[[1.8,0.45],[2.2,0.55],[2.8,0.7]]
for cp in coulombparam: 
    parms.append(
            { 
              'CONVERGED'           : 0.001,
              'FLAVORS'             : 4,
              'H'                   : 0,
              'H_INIT'              : 0.,
              'MAX_IT'              : 15,
              'MAX_TIME'            : 20,
              'MU'                  : 0,
              'N'                   : 1000,
              'NMATSUBARA'          : 1000,
              'N_MEAS'              : 10000,
              'N_ORDER'             : 50,
              'SEED'                : 0,
              'SOLVER'              : 'Hybridization',
              'SYMMETRIZATION'      : 1,
              'SWEEPS'              : 100000000,
              'BETA'                : 30,
              'THERMALIZATION'      : 10,
              'U'                   : cp[0],
              'J'                   : cp[1],
              't0'                  : 0.5,
              't1'                  : 1,
              'CHECKPOINT'          : 'dump'
        }
        )

# NOTE: this script will not be running if your ALPS repository version is older than 6243
#       ( In that case, you have to run the program from command line. )
#       WARNING: in the previous script version the parameters 't0' and 't1' were ignored (they are not supported with OMEGA_LOOP=1);
#                Thus one has simulated a totally different system!

# For more precise calculations we propose to you to:
#   enhance the MAX_TIME (to 60), 
# ( the runtime of the script with changed parameters will be roughly 45 minutes )

## Please run the tutorial5a.py before this one or uncomment the following lines.
## This tutorial relies on the results created there.

# #write the input file and run the simulation
# for p in parms:
#     input_file = pyalps.writeParameterFile('parm_u_'+str(p['U'])+'_j_'+str(p['J']),p)
#     res = pyalps.runDMFT(input_file)

listobs_ = ['Green_0', 'Green_2']   # flavor 0 is SYMMETRIZED with 1, flavor 2 is SYMMETRIZED with 3
    
ll=pyalps.load.Hdf5Loader()
for obs in listobs_:
  listobs = [obs]
  for p in parms:
    data = ll.ReadDMFTIterations(pyalps.getResultFiles(pattern='parm_u_'+str(p['U'])+'_j_'+str(p['J'])+'.h5'), measurements=listobs, verbose=True)
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
    plt.ylabel(r'$G_{flavor}(\tau)$')
    plt.title('DMFT-05: ' + r'$U = %.4s$' %nd[0].props['U'] +'; flavor='+obs[len(obs)-1])
    pyalps.plot.plot(nd)
    plt.legend()

plt.show()
