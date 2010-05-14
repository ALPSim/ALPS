import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.pyplot


#prepare the input parameters
parms=[]
for b in [6.,8.]: #,10.,12.,14.,16.
    parms.append(
            { 
              'SEED'                : 0, 
              'THERMALIZATION'      : 10000,
              'SWEEPS'              : 1000000,
              'MAX_TIME'            : 60,
              'MAX_IT'              : 1,
              'BETA'                : b,
              'SITES'               : 1,
              'N'                   : 16,
              'NMATSUBARA'          : 500, 
              'U'                   : 3,
              't'                   : 0.707106781186547,
              'MU'                  : 0,
              'H'                   : 0,
              'TOLERANCE'           : 0.0001,
              'CONVERGED'           : 0.0003,
              'SYMMETRIZATION'      : 0,
              'ANTIFERROMAGNET'     : 1,
              'SOLVER'              : '/opt/alps/bin/hirschfye',
              'FLAVORS'             : 2,
              'OMEGA_LOOP'          : 1,
              'G0OMEGA_INPUT'       : 'G0_omega_input',
              'BASENAME'            : 'hirschfye.param'
            }
        )

#write the input file and run the simulation
for p in parms:
    input_file = pyalps.writeParameterFile('parm_beta_'+str(p['BETA']),p)
    res = pyalps.runDMFT(input_file)

#data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm42a'),['Greens'])
