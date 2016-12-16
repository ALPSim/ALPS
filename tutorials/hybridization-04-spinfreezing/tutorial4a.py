from __future__ import print_function
 #############################################################################/
 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations
 #
 # ALPS Libraries
 #
 # Copyright (C) 2012 by Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 #
 #
 # This software is part of the ALPS Applications, published under the ALPS
 # Application License; you can use, redistribute it and/or modify it under
 # the terms of the license, either version 1 or (at your option) any later
 # version.
 # 
 # You should have received a copy of the ALPS Application License along with
 # the ALPS Applications; see the file LICENSE.txt. If not, the license is also
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
 #############################################################################/

 # In this tutorial, the DMFT equations are solved for a two orbital model
 # within a minimal self-consistency on the Bethe lattice.
 # The impurity model is solved using the hybridization expansion solver python
 # module. For details on how to use the solver refer to its documentation.
 #
 # Follow-up tutorials (tutorial3b.py and tutorial3c.py) illustrate how to
 # compute different quantities from the simulation results within an
 # evaluation script.
 #
 # The script can be used to reproduce the results of Figs. 4, 5 and 6 of the
 # ALPS CT-HYB paper.
 #
 # The corresponding values of the corresponding total densities and chemical
 # potentials are:
 # mu = 0.45, <n> = 0.72
 # mu = 3.60, <n> = 1.29
 # mu = 4.32, <n> = 1.53
 #
 # For simplicity, this script runs for a single set of parameters only.
 #
 # Note that this tutorial requires a significant runtime in order to obtain
 # converged results. It is recommended to run it in parallel on a cluster.
 #
 # Run this script as:
 # alpspython tutoria4a.py
 #
 # This python script is MPI aware and can hence be called using mpirun:
 #
 # mpirun -np 32 alpspython tutoria4a.py
 #
 # In case this does not work, try:
 #
 # mpirun -np 32 sh alpspython tutorial4a.py

import pyalps.mpi as mpi      # MPI library
from pyalps.hdf5 import archive   # hdf5 interface
import pyalps.cthyb as cthyb  # the solver module
from numpy import exp,sqrt,pi # some math
from numpy import array,zeros # numpy arrays
##################################################################################################################
#                                                                                                                #
#                                               P A R A M E T E R S                                              #
#                                                                                                                #
##################################################################################################################

runtime_dmft       = 600
runtime_dmft_final = 86400
dmft_iterations    = 20

# interaction parameters
U  = 8.0
J  = U/6.0
Up = U-2*J

parms = {
# solver parameters
'SWEEPS'                     : 1000000000,                         #sweeps to be done
'THERMALIZATION'             : 1000,                               #thermalization sweeps to be done
'SEED'                       : 42,                                 #random number seed
'N_MEAS'                     : 1,                                  #number of sweeps after which a measurement is done
'N_ORBITALS'                 : 4,                                  #number of 'orbitals', i.e. number of spin-orbital degrees of freedom or timelines of segments
'BASENAME'                   : "hyb.param",                        #base name of the h5 output file
'MAX_TIME'                   : runtime_dmft,                       #runtime of the solver per iteration
'VERBOSE'                    : 1,                                  #whether to output extra information
'TEXT_OUTPUT'                : 0,                                  #whether to write results in human readable (text) format
# file names
'DELTA'                      : "Delta.h5",                         #file name of the hybridization function
'DELTA_IN_HDF5'              : 1,                                  #whether to read the hybridization from an h5 archive
# physics parameters
'U'                          : U,                                  #Hubbard repulsion
"U'"                         : Up,                                 #U' parameter of the U-3J Hamiltonian
'J'                          : J,                                  #Hund's coupling J of the U-3J Hamiltonian
'MU'                         : 3.6,                                #chemical potential
'BETA'                       : 50.0,                               #inverse temperature
# measurements will be turned on in the last iterationly only
# measurement parameters
'N_HISTOGRAM_ORDERS'         : 50,                                 #maximum order for the perturbation order histogram
'N_TAU'                      : 2500,                               #number of imaginary time points (tau_0=0, tau_N_TAU=BETA)
'N_MATSUBARA'                : 512,                                #number of Matsubara frequencies
'N_nn'                       : 200,                                #number of imaginary time points for the density-density correlation function
'N_LEGENDRE'                 : 80,                                 #number of Legendre coefficients
# additional parameters (used in self-consistency only)
't'                          : 1,                                  #hopping
'mix'                        : 0.5                                 #mixing parameter for hybridization update
}# parms

if mpi.rank==0:
  print("generating initial hybridization...")
  g=[]
  I=complex(0.,1.)
  mu=0.0
  for n in range(parms['N_MATSUBARA']):
    w=(2*n+1)*pi/parms['BETA']
    g.append(2.0/(I*w+mu+I*sqrt(4*parms['t']**2-(I*w+mu)**2))) # noninteracting Green's function on Bethe lattice
  delta=[]
  for i in range(parms['N_TAU']+1):
    tau=i*parms['BETA']/parms['N_TAU']
    g0tau=0.0;
    for n in range(parms['N_MATSUBARA']):
      iw=complex(0.0,(2*n+1)*pi/parms['BETA'])
      g0tau+=((g[n]-1.0/iw)*exp(-iw*tau)).real # Fourier transform with tail subtracted
    g0tau *= 2.0/parms['BETA']
    g0tau += -1.0/2.0 # add back contribution of the tail
    delta.append(parms['t']**2*g0tau) # delta=t**2 g

  # write hybridization function to hdf5 archive (solver input)
  ar=archive(parms['DELTA'],'w')
  for m in range(parms['N_ORBITALS']):
    ar['/Delta_%i'%m]=delta
  del ar

mpi.world.barrier() #wait until hybridization is written to file

###################################################################################################################
#                                                                                                                 #
#                                       S E L F C O N S I S T E N C Y    L O O P                                  #
#                                                                                                                 #
###################################################################################################################
for it in range(dmft_iterations):

  if mpi.rank==0:
    print("****************************************************************************")
    print("*                           DMFT iteration %3i                             *"%(it))
    print("****************************************************************************")

  # !always make sure that parameters are changed on all threads equally!
  # (i.e. don't wrap this into an 'if mpi.rank==0' statement)
  if it==dmft_iterations-1:
    parms['MAX_TIME'] = runtime_dmft_final
    # turn on further measurements for final dmft interation
    parms['MEASURE_freq']=1               # Matsubara measurement (G and Sigma)
    parms['MEASURE_legendre']=1           # Legendre measurement (G and Sigma)
    parms['MEASURE_nn']=1                 # equal-time density-density correlator
    parms['MEASURE_nnt']=1                # density-density correlation function
    parms['MEASURE_sector_statistics']=1  # sector statistics
    parms['TEXT_OUTPUT']=1

  # solve the impurity model
  cthyb.solve(parms)

  # self-consistency on the master
  if mpi.rank==0:
    # read Green's function from file
    ar=archive(parms['BASENAME']+'.out.h5','r')
    # symmetrize G(tau)
    # here all orbitals and spins are degenerate
    gt=array(zeros(parms['N_TAU']+1))
    for m in range(parms['N_ORBITALS']):
      gt+=ar['G_tau/%i/mean/value'%m]
    gt/=parms['N_ORBITALS']
    del ar

    # Bethe lattice self-consistency: delta(tau)=t**2 g(tau)
    # read delta_old
    ar=archive(parms['DELTA'],'rw')
    for m in range(parms['N_ORBITALS']):
      delta_old=ar['/Delta_%i'%m]
      delta_new=array(zeros(parms['N_TAU']+1))
      delta_new=(1.-parms['mix'])*(parms['t']**2 * gt) + parms['mix']*delta_old # mix old and new delta

      # write hybridization to the h5 archive (this is solver input)
      ar['/Delta_%i'%m]=delta_new
    del ar

  mpi.world.barrier() # wait until solver input is written

# go back and loop

