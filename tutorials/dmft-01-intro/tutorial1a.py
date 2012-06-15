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
import pyalps.plot


#prepare the input parameters
parms=[]
parms.append(
            {
              'ANTIFERROMAGNET'     : 1,
              'FLAVORS'             : 2,
              'H'                   : 0.2,
              'H_INIT'              : 0.2,
              'MAX_IT'              : 3,
              'MAX_TIME'            : 60,
              'MU'                  : -0.5,
              'N'                   : 1000,
              'NMATSUBARA'          : 1000,
              'N_MEAS'              : 10000,
              'N_ORDER'             : 50,
              'OMEGA_LOOP'          : 1,
              'SEED'                : 0,
              'SITES'               : 1,
              'SOLVER'              : 'Hybridization',
              'SYMMETRIZATION'      : 0,
              'CONVERGED'           : 0.01,
              'U'                   : 2,
              't'                   : 1,
              'SWEEPS'              : 100000000,
              'THERMALIZATION'      : 1000,
              'BETA'                : 8,
              'CHECKPOINT'          : 'dump',
              'G0OMEGA_INPUT'       : 'G0_omega_noninteracting'
            }
        )
        
#write the input file and run the simulation
for p in parms:
    input_file = pyalps.writeParameterFile('parm_file',p)
    res = pyalps.runDMFT(input_file)
