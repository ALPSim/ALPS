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

 # This tutorial illustrates how to extract measured quantities from the hdf5
 # archive, perform some basic computations with them and then plot results.
 # For more information on the solver and archive paths of observables see
 # its documentation.
 #
 # This tutorial requires the results to be produced first by running the
 # tutorial4a.py script.
 #
 # Run this script as:
 # alpspython tutorial4c.py

from pyalps.hdf5 import archive # hdf5 interface
from numpy import *

ar=archive('hyb.param.out.h5','rw') # read the parameters from the output file
parms=ar['/parameters'] # returns a python dict

n_orb=parms['N_ORBITALS']

nnt = array(zeros(n_orb**2)).reshape(n_orb,n_orb)
SzSzt=[]
tau=[]
for tn in range(parms['N_nn']+1):
  tau.append(tn*parms['BETA']/parms['N_nn']) # from 0 to BETA
  for i in range(n_orb):
    for j in range(i+1): #j<=i
      nnt[i][j]=ar['/nnt_%i_%i'%(i,j)][tn]

  s=0.0
  for i in range(n_orb/2):#orb
    for j in range(n_orb/2):#orb
      for s1 in range(2):#spin
        for s2 in range(2):#spin
          idx1=2*i+s1
          idx2=2*j+s2
          sp1=1-2*s1 #s=0 -> sp=1; s=1 -> sp=-1
          sp2=1-2*s2
          if(idx2<=idx1): s+=nnt[idx1][idx2]*sp1*sp2/4.
          else: s+=nnt[idx2][idx1]*sp1*sp2/4.
  SzSzt.append(s)

ar['SzSzt']=SzSzt # save to file
del ar

# plot the results
import matplotlib.pyplot as plt
plt.figure()
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\langle S_{z}(\tau)S_z(0)\rangle$')
plt.title('hybridization-04c: Spin-spin correlation function of a two-orbital model\n(using the hybridization expansion impurity solver)')
plt.plot(tau, SzSzt)
plt.xlim(0,parms['BETA'])
plt.ylim(0,)
plt.show()



