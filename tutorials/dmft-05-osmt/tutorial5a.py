# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Brigitte Surer <surerb@phys.ethz.ch> 
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

