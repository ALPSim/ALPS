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

 # This tutorial is a minimal example illustrating the use of the python interface
 # to hybridization expansion solver.
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
 # mpirun -np 2 bash alpspython tutorial1.py

import pyalps.cthyb as cthyb # the solver module
import pyalps.mpi as mpi     # MPI library

# write a simple (constant) hybridization function to file
f=open("delta.dat","w")
for i in range(1001):
  f.write("%i %f %f\n"%(i,-0.5,-0.5))
f.close()

# specify solver parameters
parms={
'SWEEPS'              : 100000000,
'MAX_TIME'            : 10,
'THERMALIZATION'      : 1000,
'SEED'                : 0,
'N_MEAS'              : 50,
'N_HISTOGRAM_ORDERS'  : 50,
'N_ORBITALS'          : 2,
'U'                   : 5,
'MU'                  : 2.5,
'DELTA'               : "delta.dat",
'N_TAU'               : 1000,
'BETA'                : 30,
'TEXT_OUTPUT'         : 1
}

# solve the impurity model
cthyb.solve(parms)




