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

# DOS of cubic lattice

import sys
import matplotlib.pyplot as plt
from numpy import zeros
from math import sin, cos, pi, sqrt

print "Description:"
print "This program produces histogram of the density of states for cubic lattice in the tight-binding approximation, with nearest-neighbor hopping amplitude taken to be t=1."
print "by Jakub Imriska\n"
print "Enter the linear discretization GRID (cost: GRID^3/6):"
GRID = eval(raw_input('--> '))
print "  (the histogram will be computed on a grid of ", 2*GRID,' x ',2*GRID,' x ', 2*GRID,' k-points in the Brillouin zone; the program makes use of the symmetries)'
print "Enter the number of bins of the histogram (recommended divisible by 12, for Simpson integration later):"
BINS = eval(raw_input('--> '))
print "  (the output file will have ",BINS+1," rows, at the ends there are halves of bins)"

# dispersion relation is given by:
#   e(kx,ky) = -2t (\cos(k_x*a) + \cos(k_y*a) + \cos(k_z*a))
# we take t=1
# Brillouin zone: k_i * a \in <-\pi,pi)

# do histogram
lower=-6.
upper=6.
bin_width = (upper-lower) / BINS
DOS = zeros(BINS+2)

def increment(e,m):
  DOS[int(round((e - lower)/bin_width))] += m

cos_ = []
for i in range(-GRID,GRID):
  k = pi*i/GRID
  cos_.append(-2.*cos(k))

def Multiplicity(x,y,z):
  # for sorted 0<=x<=y<=z<=GRID gives the "mirror"-multiplicity
  if x>0 and z<GRID:
    return 8
  q=0
  if x==0 or x==GRID:
    q+=1
  if y==0 or y==GRID:
    q+=1
  if z==0 or z==GRID:
    q+=1
  if q==0:
    return 8
  elif q==1:
    return 4
  elif q==2:
    return 2
  else:
    return 1

def CompleteMultiplicity(x,y,z):
  # for ordinarily sorted 0<=x<=y<=z<=GRID gives the multiplicity
  if x<y and y<z:
    return 6*Multiplicity(x,y,z)
  else:
    q=0
    if x==y:
      q+=1
    if y==z:
      q+=1
    if q==1:
      return 3*Multiplicity(x,y,z)
    else:
      return Multiplicity(x,y,z)

for x in range(0,GRID+1):
  C = cos_[x]
  for y in range(x,GRID+1):
    D = C + cos_[y]
    for z in range(y,GRID+1):
      increment(D + cos_[z], CompleteMultiplicity(x,y,z))

counter=0
inc = 1./(8.*GRID*GRID*GRID*bin_width)  # normalized to 1
for x in range(0,len(DOS)):
  counter+=DOS[x]
  DOS[x]*=inc
print "Number of processed k-points: ",counter, "  (should be ", 8*GRID*GRID*GRID,')'
  
# correct normalization for the 1st and last bin
DOS[0] *= 2.      # it is a half-bin (has only half of the usual width)
DOS[BINS] *= 2.   # it is a half-bin (has only half of the usual width)

energies = zeros(BINS+2)
energies[0] = lower
energies[BINS] = upper
for j in range(0,BINS):
  energies[j] = lower + j*bin_width

def func(x,n):
  if n==0:
    return DOS[x]
  elif n==1:
    return DOS[x]*energies[x]
  else:
    return DOS[x]*energies[x]*energies[x]

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
print "  second moment of the normalized DOS = ", Integrate(2)/norm,"  (exact: 6.0)"

print "Histogram created."

plt.plot(energies[0:BINS+1],DOS[0:BINS+1],'r-')      
plt.xlabel('energy / t --->')
plt.ylabel('DOS  --->')
plt.title('DOS of the cubic lattice')
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
    file_out.write(str(DOS[j]))
    file_out.write('\n')
  file_out.close()
