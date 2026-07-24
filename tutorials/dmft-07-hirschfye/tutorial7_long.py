# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010      by Brigitte Surer <surerb@phys.ethz.ch> 
#               2012-2013 by Jakub Imriska  <jimriska@phys.ethz.ch>
# 
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.plot


#prepare the input parameters
parms=[]
for b in [6., 8., 10., 12., 14., 16.]:
    parms.append(
            { 
              'ANTIFERROMAGNET'     : 1,
              'CONVERGED'           : 0.003,
              'FLAVORS'             : 2,
              'H'                   : 0,
              'H_INIT'              : 0.05,
              'MAX_IT'              : 18,
              'MAX_TIME'            : 20,
              'MU'                  : 0,
              'N'                   : 16,
              'NMATSUBARA'          : 500, 
              'OMEGA_LOOP'          : 1,
              'SEED'                : 0, 
              'SITES'               : 1,
              'SOLVER'              : 'hirschfye',
              'SYMMETRIZATION'      : 0,
              'TOLERANCE'           : 0.001,
              'U'                   : 3,
              't'                   : 0.707106781186547,
              'SWEEPS'              : 1000000,
              'THERMALIZATION'      : 10000,
              'BETA'                : b
            }
        )

# For more precise simulation we propose to you to:
#   lower CONVERGED (to 0.0003) and TOLERANCE (to 0.0001) and increase MAX_TIME

#write the input file and run the simulation
for p in parms:
    input_file = pyalps.writeParameterFile('parm_beta_'+str(p['BETA']),p)
    res = pyalps.runDMFT(input_file)

listobs=['0', '1']
    
data = pyalps.loadMeasurements(pyalps.getResultFiles(pattern='parm_beta_*h5'), respath='/simulation/results/G_tau', what=listobs)
for d in pyalps.flatten(data):
    d.x = d.x*d.props["BETA"]/float(d.props["N"])
    d.props['label'] = r'$\beta=$'+str(d.props['BETA'])+'; flavor='+str(d.props['observable'][len(d.props['observable'])-1])
plt.figure()
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G_{flavor}(\tau)$')
plt.title('DMFT-07: Neel transition for the Hubbard model on the Bethe lattice\n(using the Hirsch-Fye impurity solver)')
pyalps.plot.plot(data)
plt.legend()
plt.show()
