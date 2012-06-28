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

#prepare the input parameters
parms=[]
parms.append(
        { 
          'ANTIFERROMAGNET'         : 0,
          'CHECKPOINT'              : 'dump_int',
          'CONVERGED'               : 0.0025,
          'CONVERGENCE_CHECK_PERIOD': 500,
          'FLAVORS'                 : 2,
          'H'                       : 0,
          'H_INIT'                  : 0.,
          'MAX_IT'                  : 12,
          'MAX_TIME'                : 120,
          'MU'                      : 0,
          'N'                       : 500,
          'NMATSUBARA'              : 500,
          'NMATSUBARA_MEASUREMENTS' : 18, 
          'NSELF'                   : 5000,
          'MEASUREMENT_PERIOD'      : 10,
          'MU'                      : 0,
          'OMEGA_LOOP'              : 1,
          'SEED'                    : 0, 
          'SITES'                   : 1,
          'SOLVER'                  : 'Interaction Expansion',
          'SYMMETRIZATION'          : 1,
          'U'                       : 3,
          't'                       : 0.707106781186547,
          'RECALC_PERIOD'           : 3000,
          'SWEEPS'                  : 100000000,
          'THERMALIZATION'          : 1000,
          'ALPHA'                   : -0.01,
          'HISTOGRAM_MEASUREMENT'   : 1,
          'BETA'                    : 32
        }
    )
    
# For more precise calculations we propose to you to:
#   enhance the MAX_TIME, MAX_IT and lower CONVERGED

#write the input file and run the simulation
input_file = pyalps.writeParameterFile('parm_int',parms[0])
res = pyalps.runDMFT(input_file)
