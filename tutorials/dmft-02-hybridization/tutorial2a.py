import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.pyplot


#prepare the input parameters
parms=[]
for b in [6.,8.,10.,12.,14.,16.]: 
    parms.append(
            { 
              'SEED'                : 0, 
              'THERMALIZATION'      : 1000,
              'SWEEPS'              : 100000000,
              'MAX_TIME'            : 60,
              'MAX_IT'              : 10,
              'BETA'                : b,
              'SITES'               : 1,
              'N'                   : 1000,
              'NMATSUBARA'          : 1000, 
              'U'                   : 3,
              't'                   : 0.707106781186547,
              'MU'                  : 0,
              'H'                   : 0,
              'H_INIT'              : 0.2,
              'TOLERANCE'           : 0.0001,
              'CONVERGED'           : 0.0003,
              'SYMMETRIZATION'      : 0,
              'ANTIFERROMAGNET'     : 1,
              'SOLVER'              : 'Hybridization',
              'FLAVORS'             : 2,
              'OMEGA_LOOP'          : 1,
              'CHECKPOINT'          : 'dump',
              'F'                   : 10,
              'N_ORDER'             : 50,
              'N_MEAS'              : 10000,
              'N_SHIFT'             : 0,
              'N_FLIP'              : 0,
              'N_MOVE'              : 0,
              'OVERLAP'             : 0
            }
        )

#write the input file and run the simulation
for p in parms:
    input_file = pyalps.writeParameterFile('parm_beta_'+str(p['BETA']),p)
    res = pyalps.runDMFT(input_file)


