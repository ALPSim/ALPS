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
 # The hybridization function is chosen such that the impurity model is equivalent
 # to a correlated site couple to a single bath site with coupling V=1 and energy
 # epsilon=0 (Delta(tau)=-V**2/2=const.) . It can therefore be compared to the
 # exact result obtained by exact diagonalization (see subdirectory ED).
 #
 # Run this script as:
 # alpspython tutorial1.py
 #

import pyalps.cthyb as cthyb # the solver module

# specify solver parameters
parms={
'SWEEPS'              : 100000000,
'MAX_TIME'            : 60,
'THERMALIZATION'      : 1000,
'SEED'                : 0,
'N_MEAS'              : 50,
'N_HISTOGRAM_ORDERS'  : 50,
'N_ORBITALS'          : 2,
'U'                   : 4.0,
'MU'                  : 2.0,
'DELTA'               : "delta.dat",
'N_TAU'               : 1000,
'BETA'                : 45,
'TEXT_OUTPUT'         : 1,
}

# Write a simple (constant) hybridization function to file.
# This corresponds to the impurity site with on-site interaction U
# coupled to a single noninteracting bath site with hybridization
# V=1 at energy epsilon=0: Delta(tau)=-V^2/2=const.

f=open("delta.dat","w")
for i in range(parms["N_TAU"]+1):
    f.write("%i %f %f\n"%(i,-0.5,-0.5))
f.close()

# solve the impurity model
cthyb.solve(parms)




