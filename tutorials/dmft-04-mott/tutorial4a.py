import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.pyplot


#prepare the input parameters
parms=[]
for u in [4.,5.,6.,8.]: 
    parms.append(
            { 
              'SEED'                    : 0, 
              'THERMALIZATION'          : 10,
              'SWEEPS'                  : 100000000,
              'MAX_TIME'                : 60,
              'MAX_IT'                  : 20,
              'BETA'                    : 20,
              'N'                       : 1000,
              'NMATSUBARA'              : 1000, 
              'NMATSUBARA_MEASUREMENTS' : 18,
              'U'                       : u,
              't'                       : 1,
              'MU'                      : 0,
              'H'                       : 0,
              'H_INIT'                  : 0.,
              'CONVERGED'               : 0.0001,
              'TOLERANCE'               : 0.0003,
              'SYMMETRIZATION'          : 0,
              'ANTIFERROMAGNET'         : 1,
              'SOLVER'                  : 'Hybridization',
              'G0OMEGA_INUT'            : 'G0_omega_input',
              'FLAVORS'                 : 2,
              'OMEGA_LOOP'              : 1,
              'F'                       : 10,
              'N_ORDER'                 : 50,
              'N_MEAS'                  : 10000,
              'N_SHIFT'                 : 0,
              'N_FLIP'                  : 0,
              'N_MOVE'                  : 0,
              'OVERLAP'                 : 0,
              'CHECKPOINT'              : 'dump'
            }
        )

#write the input file and run the simulation
for p in parms:
    input_file = pyalps.writeParameterFile('parm_u_'+str(p['U']),p)
    res = pyalps.runDMFT(input_file)
