# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2012-2013 by Jakub Imriska  <jimriska@phys.ethz.ch>
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
import matplotlib.pyplot as plt
import pyalps.plot


#prepare the input parameters
parms=[]
for u in [3.]: 
  for b in [6.]:
    parms.append(
            { 
                'BETA' : b,          # inverse temperature
                'MU' : 0.0,          # chemical potential corresponding to half-filling
                'U' : u,             # Hubbard interaction
                'FLAVORS' : 2,       # corresponds to spin up/down
                'SITES' : 1,         # number of sites in the impurity
                'H' : 0.0,           # there is no magnetic field
                'H_INIT' : 0.05,     #  we set initial field to split spin up/down in order to trigger AF phase
                'OMEGA_LOOP' : 1,        # the selfconsistency runs in Matsubara frequencies
                'ANTIFERROMAGNET' : 1,   # allow Neel order
                'SYMMETRIZATION' : 0,    # do not enforce paramagnetic solution
                'NMATSUBARA' : 500,      # number of Matsubara frequencies
                'N' : 500,               # bins in imaginary time
                'CONVERGED' : 0.005,     # criterion for convergency
                'MAX_TIME' : 60,         # max. time spent in solver in a single iteration in seconds
                'G0OMEGA_INPUT' : "",    # forces to start from the local non-interacting Green's function
                'MAX_IT' : 10,           # max. number of self-consistency iterations
                'SWEEPS' : 10000,    # max. number of sweeps in a single iteration
                'THERMALIZATION' : 500, # number of thermalization sweeps
                'SEED' : 0,              # random seed
                'SOLVER' : "hybridization",   # we take the hybridization impurity solver
                'SC_WRITE_DELTA' : 1,         # input for the hybridization impurity solver is the hybridization function Delta, which has to be written by the selfconsistency
                'N_MEAS' : 5000,              # number of Monte Carlo steps between measurements
                'N_ORDER' : 50,               # histogram size
                'TWODBS' : 1,     # the Hilbert transformation integral runs in k-space, sets square lattice
                't' : 1,          # the nearest-neighbor hopping
                'tprime' : 0,     # the second nearest-neighbor hopping
                'L' : 64,         # discretization in k-space in the Hilbert transformation
                'GENERAL_FOURIER_TRANSFORMER' : 1,  # Fourier transformer for a general bandstructure
                'EPS_0' : 0,                        # potential shift for the flavor 0
                'EPS_1' : 0,                        # potential shift for the flavor 1
                'EPSSQ_0' : 4,                      # the second moment of the bandstructure for the flavor 0
                'EPSSQ_1' : 4,                      # the second moment of the bandstructure for the flavor 1
            }
        )

#write the input file and run the simulation
for p in parms:
    input_file = pyalps.writeParameterFile('hybrid_TWODBS_beta_'+str(p['BETA'])+'_U_'+str(p['U']),p)
    res = pyalps.runDMFT(input_file)

listobs=['0']  # we look only at flavor=0
    
data = pyalps.loadMeasurements(pyalps.getResultFiles(pattern='hybrid_TWODBS*h5'), respath='/simulation/results/G_tau', what=listobs, verbose=True)
for d in pyalps.flatten(data):
    d.x = d.x*d.props["BETA"]/float(d.props["N"])
    d.props['label'] = r'$\beta=$'+str(d.props['BETA'])
plt.figure()

plt.xlabel(r'$\tau$')
plt.ylabel(r'$G_{flavor=0}(\tau)$')
plt.title('DMFT-08, TWODBS option: Hubbard model on the square lattice')
pyalps.plot.plot(data)
plt.legend()
plt.show()
