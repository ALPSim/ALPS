import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.pyplot


#prepare the input parameters
parms=[]
for b in [6.,8.,10.,12.,14.,16.]: 
    parms.append(
            { 
              'SEED'                    : 0, 
              'THERMALIZATION'          : 1000,
              'SWEEPS'                  : 100000,
              'NRUNS'                   : 1,
              'MAX_TIME'                : 60,
              'MAX_IT'                  : 10,
              'MEASUREMENT_PERIOD'      : 10,
              'RECALC_PERIOD'           : 300,
              'CONVERGENCE_CHECK_PERIOD': 500,
              'BETA'                    : b,
              'SITES'                   : 1,
              'N'                       : 500,
              'NMATSUBARA'              : 500, 
              'NMATSUBARA_MEASUREMENTS' : 18,
              'U'                       : 3,
              't'                       : 0.707106781186547,
              'MU'                      : 0,
              'H'                       : 0,
              'H_INIT'                  : 0.,
              'CONVERGED'               : 0.0008,
              'SYMMETRIZATION'          : 0,
              'ANTIFERROMAGNET'         : 1,
              'SOLVER'                  : 'Interaction Expansion',
              'HISTOGRAM_MEASUREMENT'   : 1,
              'NSELF'                   : 5000,
              'ALPHA'                   : -0.01,
              'G0OMEGA_INUT'            : 'G0_omega_input',
              'FLAVORS'                 : 2,
              'OMEGA_LOOP'              : 1,
              'CHECKPOINT'              : 'dump',
            }
        )

#write the input file and run the simulation
for p in parms:
    input_file = pyalps.writeParameterFile('parm_beta_'+str(p['BETA']),p)
    res = pyalps.runDMFT(input_file)
