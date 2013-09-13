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
for b in [6.,8.,10.,12.,14.,16.]:
    parms.append(
            {
              'ANTIFERROMAGNET'     : 1,
              'CONVERGED'           : 0.005,
              'FLAVORS'             : 2,
              'H'                   : 0,
              'H_INIT'              : 0.05,
              'MAX_IT'              : 15,
              'MAX_TIME'            : 15,
              'MU'                  : 0,
              'N'                   : 500,
              'NMATSUBARA'          : 500,
              'N_MEAS'              : 10000,
              'N_ORDER'             : 50,
              'OMEGA_LOOP'          : 1,
              'SEED'                : 0,
              'SITES'               : 1,
              'SOLVER'              : 'hybridization',  # name of the solver executable, evt. full path needed
              'SC_WRITE_DELTA'      : 1,
              'SYMMETRIZATION'      : 0,
              'U'                   : 3,
              't'                   : 0.707106781186547,
              'SWEEPS'              : 100000000,
              'THERMALIZATION'      : 1000,
              'BETA'                : b,
              'CHECKPOINT'          : 'solverdump_beta_'+str(b)
            }
        )

# NOTE: in revision of ALPS older than 6238, the MAX_TIME will effectively be 60 seconds.        
# For more precise calculations we propose to you to:
#   enhance the MAX_TIME (to 60), 
#   lower the CONVERGED (to 0.003), 
#   increase MAX_IT (to 20)
#   raise N and NMATSUBARA (to 1000)
# ( the runtime of the script with changed parameters will be roughly 2 hours )
        
#write the input file and run the simulation
for p in parms:
    input_file = pyalps.writeParameterFile('parm_beta_'+str(p['BETA']),p)
    res = pyalps.runDMFT(input_file)

flavors=parms[0]['FLAVORS']
listobs=[]   
for f in range(0,flavors):
    listobs.append(str(f))
    
data = pyalps.loadMeasurements(pyalps.getResultFiles(pattern='parm_beta_*h5'), respath='/simulation/results/G_tau', what=listobs, verbose=True)
for d in pyalps.flatten(data):
    d.x = d.x*d.props["BETA"]/float(d.props["N"])
    d.props['label'] = r'$\beta=$'+str(d.props['BETA'])+'; flavor='+str(d.props['observable'][len(d.props['observable'])-1])

plt.figure()
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G_{flavor}(\tau)$')
plt.title('DMFT-02: Neel transition for the Hubbard model on the Bethe lattice\n(using the Hybridization expansion impurity solver)')
pyalps.plot.plot(data)
plt.legend()
plt.show()
