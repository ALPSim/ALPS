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

 # This tutorial implements the DMFT self-consistency on the Bethe-lattice
 # for the Hubbard-Holstein model with a single phonon mode at w0 and coupling
 # lambda. The latter gives rise to a retarded interaction in the impurity model.
 # The script can be used to reproduce the results of Fig.3 of the ALPS CT-HYB paper.
 #
 # The impurity model is solved using the cthyb module (python interface to the
 # hybridization solver). For details on how to use the solver see its documentation.
 #
 # The tutorial shows how to flexibly implement a selfconsistency loop without using
 # the ALPS DMFT framework. The selfconsistency is based on G(tau) and additional
 # quantities (such as self-energy etc.) are measured in the last iteration only.
 # These are written to stdout in human-readable (plottable) format.
 #
 # Run this script as:
 # alpspython tutorial2.py
 #
 # This python script is MPI aware and can hence be called using mpirun:
 #
 # mpirun -np 2 alpspython tutorial2.py
 #
 # In case this does not work, try:
 #
 # mpirun -np 2 sh alpspython tutorial2.py

import shutil
import pyalps.mpi as mpi                # mpi library
from pyalps.hdf5 import archive            # hdf5 interface
import pyalps.cthyb as cthyb            # the solver module
from numpy import sqrt,cosh,sinh,exp,pi #some math
from numpy import array,zeros,append
##################################################################################################################
#                                                                                                                #
#                                               P A R A M E T E R S                                              #
#                                                                                                                #
##################################################################################################################

# This script takes some time to run to get converged results. Make sure the runtime is large enough to get
# sensible results (depends on the number of processes).
runtime_dmft=60        # runtime for each DMFT iteration
runtime_dmft_final=600 # increase runtime on final iteration for additional measurements
dmft_iterations=10     # number of DMFT iterations

# perform calculations for fixed U and Uscr and given w0
U      = 8.0
Uscr   = 3.0
w0     = 12.0
Lambda = sqrt((U-Uscr)*w0/2.0) # choose lambda as to get the given Uscr = U - 2*lambda^2/w0

# for simplicity, this script is for a single set of parameters only
# we list all solver parameters here for completeness
parms = {
# solver parameters
# general
'SWEEPS'                     : 1000000000,                         #sweeps to be done
'THERMALIZATION'             : 1000,                               #thermalization sweeps to be done
'SEED'                       : 42,                                 #random number seed
'N_MEAS'                     : 100,                                #number of sweeps after which a measurement is done
'N_ORBITALS'                 : 2,                                  #number of 'orbitals', i.e. number of spin-orbital degrees of freedom or segments
'BASENAME'                   : "hyb.param",                        #base name of the h5 output file
'MAX_TIME'                   : runtime_dmft,                       #runtime of the solver per iteration
'VERBOSE'                    : 1,                                  #whether to output extra information
'COMPUTE_VERTEX'             : 0,                                  #whether to compute the vertex function
'TEXT_OUTPUT'                : 0,                                  #whether to write results in human readable (text) format
# file names
'DELTA'                      : "Delta.h5",                         #file name of the hybridization function
'DELTA_IN_HDF5'              : 1,                                  #whether to read the hybridization from an h5 archive
'RET_INT_K'                  : "K_tau.h5",                         #file name of retarded interaction function
'K_IN_HDF5'                  : 1,                                  #whether to read the retarded interaction from an h5 archive
# physics parameters
'U'                          : U,                                  #Hubbard repulsion
'MU'                         : U/2.-2*Lambda**2/w0,                #chemical potential (MU=U/2-2K'(0) corresponds to half-filling; here K(0)=Lambda^2/w0)
'BETA'                       : 50.0,                               #inverse temperature
# measurements
'MEASURE_freq'               : 0,                                  #whether to measure single-particle Green's function on Matsubara frequencies
'MEASURE_legendre'           : 0,                                  #whether to measure single-particle Green's function in Legendre polynomial basis
'MEASURE_g2w'                : 0,                                  #whether to measure two-particle Green's function on Matsubara frequencies
'MEASURE_h2w'                : 0,                                  #whether to measure the higher-order correlation function for the vertex on Matsubara frequencies
'MEASURE_nn'                 : 0,                                  #whether to measure equal-time density-density correlations
'MEASURE_nnt'                : 0,                                  #whether to measure the density-density correlation function (local susceptibility) in imaginary time
'MEASURE_nnw'                : 0,                                  #whether to measure the density-density correlation function (local susceptibility) on Matsubara frequencies
'MEASURE_sector_statistics'  : 0,                                  #whether to measure sector statistics
# measurement parameters
'N_HISTOGRAM_ORDERS'         : 50,                                 #maximum order for the perturbation order histogram
'N_TAU'                      : 5000,                               #number of imaginary time points (tau_0=0, tau_N_TAU=BETA)
'N_MATSUBARA'                : 512,                                #number of Matsubara frequencies
'N_nn'                       : 5000,                               #number of imaginary time points for the density-density correlation function
'N_W'                        : 20,                                 #number of bosonic Matsubara frequencies for the two-particle Green's function or local susceptibility
'N_w2'                       : 20,                                 #number of fermionic Matsubara frequencies for the two-particle Green's function
'N_LEGENDRE'                 : 80,                                 #number of Legendre coefficients
# additional parameters (used outside the solver only)
't'                          : 1,                                  #hopping
'Uscr'                       : Uscr,                               #screened interaction
'lambda'                     : Lambda,                             #electron phonon coupling
'w0'                         : w0,                                 #phonon frequency
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

  print("generating retarded interaction...")
  l   =parms['lambda']
  w0  =parms['w0']
  beta=parms['BETA']

  K  = lambda tau: - (l**2)*(cosh(w0*(beta/2.0-tau))/sinh(w0*beta/2.0) - cosh(w0*beta/2.0)/sinh(w0*beta/2.0) )/(w0*w0)
  Kp = lambda tau: + (l**2)*(sinh(w0*(beta/2.0-tau))/sinh(w0*beta/2.0))/w0

  print("U_screened =", parms['U'] - 2*l**2/w0)

  k_tau=[]
  kp_tau=[]
  for i in range(parms['N_TAU']+1):
    tau=i*parms['BETA']/parms['N_TAU']
    k_tau.append(K(tau))
    kp_tau.append(Kp(tau))

  # write retarded interaction function K(tau) and its derivative to file (solver input)
  ar=archive(parms['RET_INT_K'],'w')
  ar['/Ret_int_K']=k_tau
  ar['/Ret_int_Kp']=kp_tau
  del ar

  if(parms['TEXT_OUTPUT']==1):
    f=open('Ktau.dat','w')
    for i in range(len(k_tau)):
      tau=i*parms['BETA']/parms['N_TAU']
      f.write("%f %f %f\n"%(tau,k_tau[i],kp_tau[i]))

mpi.world.barrier() # wait until solver input is written to file

###################################################################################################################
#                                                                                                                 #
#                                D M F T   S E L F C O N S I S T E N C Y    L O O P                               #
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
    # turn on additional measurements for the final dmft interation
    parms['MEASURE_freq']=1 # turn of Matsubara measurement
    parms['MEASURE_legendre']=1 # turn on Legendre measurement
    parms['TEXT_OUTPUT']=1  # this will write results of the final iteration in text format

  # write parameters for reference (on master only)
  if mpi.rank==0:
    ar=archive(parms['BASENAME']+'.h5','a')
    ar['/parameters']=parms
    ar['/parameters%i'%it]=parms # this is a backup for each iteration
    del ar

  # solve the impurity model in parallel
  cthyb.solve(parms)

  # self-consistency on the master
  if mpi.rank==0:
    if(parms['TEXT_OUTPUT']==1):
      shutil.copy("Gt.dat", "Gt%i.dat"%it) # keep Green's function for each iteration for monitoring
      shutil.copy("simulation.dat", "simulation%i.dat"%it) # keep some basic information for each iteration

    # read Green's function from file
    ar=archive(parms['BASENAME']+'.out.h5','r')
    # symmetrize G(tau)
    # in this case all orbitals and spins are degenerate
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

    # write input hybridization function for reference
    if(parms['TEXT_OUTPUT']==1):
      f=open('Delta%i.dat'%it,'w')
      for i in range(parms['N_TAU']+1):
        f.write("%f"%(i*parms['BETA']/parms['N_TAU']))
        for m in range(parms['N_ORBITALS']):
          f.write(" %f"%delta_old[i])
        f.write("\n")
      f.close()

  mpi.world.barrier() # wait until solver input is written

# go back and loop

