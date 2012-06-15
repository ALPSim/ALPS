# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2012 by Jakub Imriska <jimriska@phys.ethz.ch> 
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

# DOS of Hexagonal lattice

import sys
import matplotlib.pyplot as plt
from numpy import zeros
from math import sin, cos, pi, sqrt

print "Description:"
print "This program produces histogram of the density of states for hexagonal lattice in the tight-binding approximation, with hopping amplitude taken to be t=1."
print "by Jakub Imriska\n"
print "Enter the linear discretization GRID (cost: GRID^2):"
GRID = eval(raw_input('--> '))
print "  (the histogram will be computed on a grid of ", 2*GRID,' x ',2*GRID,' k-points in the Brillouin zone; the program makes use of the symmetries)'
print "Enter the number of bins of the histogram (recommended divisible by 12, for Simpson integration later):"
BINS = eval(raw_input('--> '))
print "  (the output file will have ",BINS+1," rows, at the ends there are halves of bins)"

# According to notes, the BZ may be chosen as:
#    k_x \in <0, 4\pi/3);  k_y \in <0, 2\pi/\sqrt{3})
# and the dispersion is then 
#    epsilon = \pm t \sqrt{1+4\cos\tilde{k_y}(\cos\tilde{k_y} + \cos\tilde{k_x})}
# with 
#    \tilde{k_x} = 3/2*k_x;  \tilde{k_y} = \sqrt{3}/2*k_y
# Thus, we work here directly in the \tilde{k_i} coordinates and the BZ is chosen to be symmetric
#    \tilde{k_x} \in <-\pi,\pi);  \tilde{k_y} \in <-\pi/2,\pi/2)
#
# The hopping amplitude is taken to be
#   t=1

cos_x = zeros(GRID+2)
cos_y = zeros(GRID+2)
for i in range(0,GRID+1):
  cos_x[i] = cos(pi*i/GRID)
  cos_y[i] = cos(pi*(i/(2.*GRID)))

print " ... histogram building in progres ..."

# do histogram
lower=0.  # we look only on the plus branch and the minus branch will be added at the end
upper=3. 
bin_width = (upper-lower) / (BINS/2)
DOS = zeros(BINS/2+2)

def increment(kx,ky,m):
  e = sqrt(1. + 4.*cos_y[ky]*(cos_y[ky]+cos_x[kx]))
  DOS[int(round((e - lower)/bin_width))] += m

increment(0,0,1)
increment(GRID,0,1)
increment(0,GRID,1)
increment(GRID,GRID,1)

for kx in range(1,GRID):
  increment(kx,0,2)
  increment(kx,GRID,2)

for ky in range(1,GRID):
  increment(0,ky,2)
  increment(GRID,ky,2)

for kx in range(1,GRID):
  for ky in range(1,GRID):
    increment(kx,ky,4)


counter=0
inc = 1./(4.*GRID*GRID*bin_width)/2.  # normalized to 1, both bands
for x in range(0,len(DOS)):
  counter+=DOS[x]
  DOS[x]*=inc
print "Number of processed k-points: ",counter, "  (should be ", 4*GRID*GRID,')'

  
# correct normalization for the 1st and last bin
DOS[0] = 0    # the precise value; do not set to zero, if you do not want to use Simpson integration over the DOS
DOS[BINS/2] *= 2.  # it is a half-bin (has only half of the usual width)

def get_dos(x):
  if x<0:
    return DOS[-x]
  else:
    return DOS[x]

energies, dos = zeros(BINS+2), zeros(BINS+2)
for j in range(0,BINS+1):
  energies[j] = -upper + j*bin_width
  dos[j] = get_dos(j-BINS/2)
energies[BINS] = upper

def func(x,n):
  if n==0:
    return dos[x]
  elif n==1:
    return dos[x]*energies[x]
  else:
    return dos[x]*energies[x]*energies[x]

def Integrate(n):
  if ( BINS % 2 == 1 ):
    return " error: for Simpson integration use even number of bins"
  sum1 = 0
  sum2 = 0
  halfstep = bin_width
  
  for i in range(1, BINS-1, 2):
    sum2 += func(i,n);
    sum1 += func(i+1,n);
  
  sum1 = 2. * sum1 + 4. * (sum2+func(BINS-1,n)) + func(0,n) + func(BINS,n);
  return sum1 * halfstep / 3.

print "Checks:"
norm = Integrate(0)
print "  normalization = ", norm,"  (close to 1)"
print "  first moment of the normalized DOS = ", Integrate(1)/norm,"  (exact: 0.0)"
print "  second moment of the normalized DOS = ", Integrate(2)/norm,"  (exact: 3.0)"

print "Histogram created."

plt.plot(energies[0:BINS+1],dos[0:BINS+1],'r-')      
plt.xlabel('energy / t --->')
plt.ylabel('DOS  --->')
plt.title('DOS of the hexagonal lattice')
plt.show()

print "Do you wish to save the histogram [y/n] ?"
answer = raw_input('--> ')

if answer[0]=='y':
  # write into file
  print "Set the name for the histogram output file:"
  file_name = raw_input('--> ')
  file_out = open(file_name,'w')
  for j in range(0, BINS+1):
    file_out.write(str(energies[j]))
    file_out.write('  ')
    file_out.write(str(dos[j]))
    file_out.write('\n')
  file_out.close()
