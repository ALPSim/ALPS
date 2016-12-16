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

 # This tutorial illustrates how to extract measured quantities from the hdf5
 # archive, perform some basic computations with them and then plot results.
 #
 # This tutorial requires the results to be produced first by running the
 # tutorial4a.py script.
 #
 # Run this script as:
 # alpspython tutoria4b.py

from pyalps.hdf5 import archive #hdf5 interface
from numpy import *

# function to compute the U matrix from the interaction parameters
def generate_u_matrix(FLAVORS, U, Up, J):
  udata=array(zeros(FLAVORS*FLAVORS))
  for i in range(0,FLAVORS,2):
    udata[ i   *FLAVORS+ i   ]=0 # Pauli
    udata[(i+1)*FLAVORS+(i+1)]=0 # Pauli
    udata[ i   *FLAVORS+(i+1)]=U # Hubbard repulsion same band
    udata[(i+1)*FLAVORS+ i   ]=U # Hubbard repulsion same band
    for j in range(0,FLAVORS,2):
      if(j!=i):
        udata[ i   *FLAVORS+ j   ]=Up-J # Hubbard repulsion interband same spin
        udata[(i+1)*FLAVORS+(j+1)]=Up-J # Hubbard repulsion interband same spin
        udata[ i   *FLAVORS+(j+1)]=Up # Hubbard repulsion interband opposite spin
        udata[(i+1)*FLAVORS+ j   ]=Up # Hubbard repulsion interband opposite spin
  return udata.reshape(FLAVORS,FLAVORS)

ar=archive('hyb.param.out.h5','rw') #read the parameters from the output file
parms=ar['/parameters'] # returns a python dict

n_orb = parms['N_ORBITALS']

a=array(zeros(n_orb**2)).reshape(n_orb,n_orb) # density-density correlation matrix
for i in range(1,n_orb): # i>j
  for j in range(i):
    a[i][j]=ar['/simulation/results/nn_%i_%i/mean/value'%(i,j)]

for i in range(n_orb): # i==j
  a[i][i]=ar['/simulation/results/density_%i/mean/value'%i]

for j in range(1,n_orb): # the matrix is symmetric
  for i in range(j): # i<j
    a[i][j]=a[j][i]

u=generate_u_matrix(parms['N_ORBITALS'],parms['U'],parms["U'"],parms['J'])

print('U matrix:')
print(u)

print('density-density correlation matrix:')
print(a)


# calculate moments of the self-energy high-frequency tail from density-desnsity correlators
# Sigma_tail(i\omega_n) ~ sigma0 + sigma1/(i*omega_n)
sigma0 = 0.
sigma1 = 0.

for i in range(n_orb): # average over all orbitals i
  for j in range(n_orb):
    sigma0 += u[i][j]*a[j][j]
sigma0/=n_orb
print("sigma0=", sigma0) # print average

for i in range(n_orb): # average over all orbitals i
  for k in range(n_orb):
    for l in range(n_orb):
      sigma1 += u[i][k]*u[i][l]*(a[k][l]-a[k][k]*a[l][l])
sigma1/=n_orb
print("sigma1=", sigma1) # print average

N=0.
for i in range(n_orb):
  N+=ar['/simulation/results/density_%i/mean/value'%(i)]

print("total density=" , N)

Xvalues=[]
Tvalues=[]

for n in range(parms['N_MATSUBARA']):
  Iw=complex(0.,(2*n+1)*pi/parms['BETA'])
  Xvalues.append(Iw.imag)
  Tvalues.append((sigma0+sigma1/Iw).imag) # self-energy tail

Sw=array(zeros(parms['N_MATSUBARA']), dtype=complex)
Swl=array(zeros(parms['N_MATSUBARA']), dtype=complex)

for m in range(parms['N_ORBITALS']):
  Sw+=ar['S_omega/%i/mean/value'%m]/parms['N_ORBITALS']
  Swl+=ar['S_l_omega/%i/mean/value'%m]/parms['N_ORBITALS']

del ar # close the archive

# plot the results
import matplotlib.pyplot as plt
plt.figure()
plt.xlabel(r'$\omega_n$')
plt.ylabel(r'Im$\Sigma_{ii}(i\omega_n)$')
plt.title('hybridization-04b: Self-energy of a two-orbital model\n(using the hybridization expansion impurity solver)')
plt.plot(Xvalues, imag(Sw),  label='Im $\Sigma(i\omega_n)$ (Mats.)')
plt.plot(Xvalues, imag(Swl), label='Im $\Sigma(i\omega_n)$ (Leg.)')
plt.plot(Xvalues, Tvalues, label="Tail")
plt.xlim(0,60)
plt.ylim(-2,0)
plt.legend(loc='lower right')
plt.show()





















