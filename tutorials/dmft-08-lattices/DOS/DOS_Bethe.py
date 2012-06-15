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

# DOS of the Bethe lattice in infinitely many dimensions

import sys
import matplotlib.pyplot as plt
from numpy import zeros
from math import pi, sqrt

print "Description:"
print "This program produces file containing the density of states for the Bethe lattice (in infinitely many dimensions) to be read as input by ALPS application DMFT (with option DOSFILE)."
print "by Jakub Imriska\n"
print "Enter the discretization of the density of states (recommended divisible by 2, for Simpson integration later):"
BINS = eval(raw_input('--> '))
print "  (the output file will have ",BINS+1," rows, at the ends there are halves of bins)"

# The DOS of Bethe lattice is semicircular, thus the normalized DOS is simply
#   DOS(e) = 
# we take t=1; D=2W=4t

DOS_WIDTH = 2.
NORMALIZATION = 1./(2.*pi)

energies = zeros(BINS+2)
dos = zeros(BINS+2)
dos[0] = dos[BINS] = 0.
energies[0] = -DOS_WIDTH
energies[BINS] = DOS_WIDTH

bin_width = 2.*DOS_WIDTH / BINS

for i in range(1,BINS):
  energies[i] = energies[0] + bin_width * i
  dos[i] = NORMALIZATION * sqrt(DOS_WIDTH*DOS_WIDTH - energies[i]*energies[i])

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
print "  second moment of the normalized DOS = ", Integrate(2)/norm,"  (exact: 1.0)"

plt.plot(energies[0:BINS+1],dos[0:BINS+1],'r-')      
plt.xlabel('energy / t --->')
plt.ylabel('DOS  --->')
plt.title('DOS of the Bethe lattice')
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
