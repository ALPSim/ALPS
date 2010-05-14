import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.pyplot


#prepare the input parameters
parms=[]
coulombparam=[[1.8,0.45],[2.2,0.55],[2.8,0.7]]
for cp in coulombparam: 
    parms.append(
            { 
              'SWEEPS'              : 100000000,
              'THERMALIZATION'      : 10,
              'MAX_IT'              : 20,
              'MAX_TIME'            : 180,
              'BETA'                : 30,
              'MU'                  : 0,
              'H'                   : 0,
              'H_INIT'              : 0,
              'U'                   : cp[0],
              'J'                   : cp[1],
              't'                   : 1,
              't0'                  : 0.5,
              't1'                  : 1,
              'SYMMETRIZATION'      : 1,
              'N'                   : 1000,
              'NMATSUBARA'          : 1000,
              'FLAVORS'             : 4,
              'CONVERGED'           : 0.0001,
              'TOLERANCE'           : 0.003,
              'SOLVER'              : 'Hybridization',
              'SEED'                : 0,
              'F'                   : 10,
              'N_ORDER'             : 50,
              'N_MEAS'              : 10000,
              'N_SHIFT'             : 0,
              'N_FLIP'              : 0,
              'N_MOVE'              : 0,
              'OVERLAP'             : 0,
              'CHECKPOINT'          : 'dump',
              'G0TAU_INPUT'         :'G0_tau_input'
            }
        )

#write the input file and run the simulation
for p in parms:
    input_file = pyalps.writeParameterFile('parm_u_'+str(p['U'])+'_j_'+str(p['J']),p)
    res = pyalps.runDMFT(input_file)

