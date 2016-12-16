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

 # This tutorial shows how to repeatedly solve a single impurity Anderson model
 # using the hybridization expansion solver and how to extract the calculated
 # observables and their error via the python interface.
 #
 # The script can be used to reproduce the results of Figs. 1 and 2 of the ALPS
 # CT-HYB paper.
 #
 # The results show the decrease of the effective local moment of the impurity
 # with decreasing temperature due to Kondo screening. For simplicity, a hybri-
 # dization function corresponding to a semielliptical density of states is used.
 #
 # The impurity model is solved using the cthyb module (python interface to the
 # hybridization solver). For details on how to use the solver see the documentation
 #
 # Run this script as:
 # alpspython tutorial1.py
 #
 # This python script is MPI aware and can hence be called using mpirun:
 #
 # mpirun -np 2 alpspython tutorial1.py
 #
 # In case this does not work, try:
 #
 # mpirun -np 2 sh alpspython tutorial1.py

import pyalps.mpi as mpi          # MPI library
from pyalps.hdf5 import archive       # hdf5 interface
import pyalps.cthyb as cthyb      # the solver module
import matplotlib.pyplot as plt   # for plotting results
from numpy import exp,log,sqrt,pi # some math

# generate a sequence of temperatures between a and b which are equidistant on a logarithmic scale
N_T  = 10    # number of temperatures
Tmin = 0.05  # maximum temperature
Tmax = 100.0 # minimum temperature
Tdiv = exp(log(Tmax/Tmin)/N_T)
T=Tmax
Tvalues=[]
for i in range(N_T+1):
  Tvalues.append(T)
  T/=Tdiv

Uvalues=[0.,2.] # the values of the on-site interaction
N_TAU = 1000    # number of tau-points; must be large enough for the lowest temperature (set to at least 5*BETA*U)
runtime = 1     # solver runtime (in seconds)

values=[[] for u in Uvalues]
errors=[[] for u in Uvalues]
for un,u in enumerate(Uvalues):
  parameters=[]
  for t in Tvalues:
      # prepare the input parameters; they can be used inside the script and are passed to the solver
      parameters.append(
             {
               # solver parameters
               'SWEEPS'                     : 1000000000,                         # sweeps to be done
               'THERMALIZATION'             : 1000,                               # thermalization sweeps to be done
               'SEED'                       : 42,                                 # random number seed
               'N_MEAS'                     : 10,                                 # number of sweeps after which a measurement is done
               'N_ORBITALS'                 : 2,                                  # number of 'orbitals', i.e. number of spin-orbital degrees of freedom or segments
               'BASENAME'                   : "hyb.param_U%.1f_BETA%.3f"%(u,1/t), # base name of the h5 output file
               'MAX_TIME'                   : runtime,                            # runtime of the solver per iteration
               'VERBOSE'                    : 1,                                  # whether to output extra information
               'TEXT_OUTPUT'                : 0,                                  # whether to write results in human readable (text) format
               # file names
               'DELTA'                      : "Delta.h5",                         # file name of the hybridization function
               'DELTA_IN_HDF5'              : 1,                                  # whether to read the hybridization from an h5 archive
               # physical parameters
               'U'                          : u,                                  # Hubbard repulsion
               'MU'                         : u/2.,                               # chemical potential
               'BETA'                       : 1/t,                                # inverse temperature
               # measurements
               'MEASURE_nnw'                : 1,                                  # measure the density-density correlation function (local susceptibility) on Matsubara frequencies
               'MEASURE_time'               : 0,                                  # turn of imaginary-time measurement
               # measurement parameters
               'N_HISTOGRAM_ORDERS'         : 50,                                 # maximum order for the perturbation order histogram
               'N_TAU'                      : N_TAU,                              # number of imaginary time points (tau_0=0, tau_N_TAU=BETA)
               'N_MATSUBARA'                : int(N_TAU/(2*pi)),                  # number of Matsubara frequencies
               'N_W'                        : 1,                                  # number of bosonic Matsubara frequencies for the local susceptibility
               # additional parameters (used outside the solver only)
               't'                          : 1,                                  # hopping
             }
          )

  for parms in parameters:

    if mpi.rank==0:
      ar=archive(parms['BASENAME']+'.out.h5','a')
      ar['/parameters']=parms
      del ar
      print("creating initial hybridization...") 
      g=[]
      I=complex(0.,1.)
      mu=0.0
      for n in range(parms['N_MATSUBARA']):
        w=(2*n+1)*pi/parms['BETA']
        g.append(2.0/(I*w+mu+I*sqrt(4*parms['t']**2-(I*w+mu)**2))) # use GF with semielliptical DOS
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

    mpi.world.barrier() # wait until hybridization is written to file

    # solve the impurity model in parallel
    cthyb.solve(parms)

    if mpi.rank==0:
      # extract the local spin susceptiblity
      ar=archive(parms['BASENAME']+'.out.h5','w')
      nn_0_0=ar['simulation/results/nnw_re_0_0/mean/value']
      nn_1_1=ar['simulation/results/nnw_re_1_1/mean/value']
      nn_1_0=ar['simulation/results/nnw_re_1_0/mean/value']
      dnn_0_0=ar['simulation/results/nnw_re_0_0/mean/error']
      dnn_1_1=ar['simulation/results/nnw_re_1_1/mean/error']
      dnn_1_0=ar['simulation/results/nnw_re_1_0/mean/error']

      nn  = nn_0_0 + nn_1_1 - 2*nn_1_0
      dnn = sqrt(dnn_0_0**2 + dnn_1_1**2 + ((2*dnn_1_0)**2) )

      ar['chi']=nn/4.
      ar['dchi']=dnn/4.

      del ar
      T=1/parms['BETA']
      values[un].append(T*nn[0])
      errors[un].append(T*dnn[0])

# go back and loop

if mpi.rank==0: #allow other nodes to exit
  plt.figure()
  plt.xlabel(r'$T$')
  plt.ylabel(r'$4T\chi_{dd}$')
  plt.title('hybridization-02: Kondo screening of an impurity\n(using the hybridization expansion impurity solver)')
  for un in range(len(Uvalues)):
    plt.errorbar(Tvalues, values[un], yerr=errors[un], label="U=%.1f"%Uvalues[un])
  plt.xscale('log')
  plt.legend()
  plt.show()




