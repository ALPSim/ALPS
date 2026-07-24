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
              'H_INIT'              : 0.03*b/8.,
              'MAX_IT'              : 6,
              'MAX_TIME'            : 300,
              'MU'                  : 0,
              'N'                   : 250,
              'NMATSUBARA'          : 250,
              'N_MEAS'              : 10000,
              'OMEGA_LOOP'          : 1,
              'SEED'                : 0,
              'SITES'               : 1,
              'SOLVER'              : 'hybridization',
              'SC_WRITE_DELTA'      : 1,
              'SYMMETRIZATION'      : 0,
              'U'                   : 3,
              't'                   : 0.707106781186547,
              'SWEEPS'              : int(10000*b/16.),
              'THERMALIZATION'      : 1000,
              'BETA'                : b
            }
        )


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
plt.title('DMFT-02: Neel transition for the Hubbard model on the Bethe lattice\n(using the Hybridization expansion impurity solver)')
pyalps.plot.plot(data)
plt.legend()
plt.show()
