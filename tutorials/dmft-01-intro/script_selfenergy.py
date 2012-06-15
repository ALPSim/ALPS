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
#
# WARNING: Does work with python library scripts of revision 6195 and higher. 
#


import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.plot
from math import pi

print "DESCRIPTION: This script shows the iteration-resolved selfenergy from the specifyied result file. The selfenergy is calculated via Dyson equation from the stored G and G0."
print

print "Enter the name (prefix) of the (single) result file, which is to examine :"
res_file = raw_input('--> ')

tau_repr = False

real_part = False
print "Do you want to plot the imaginary part of the selfenergy [default] or the real part [type: 'r'] ?"
answer2 = raw_input('--> ')
if len(answer2)>0 and answer2[0]=='r':
  real_part = True

print "Please select flavor [default: 0] :"
answer4 = raw_input('--> ')
flavor = 0
if len(answer4)>0:
  flavor = eval(answer4)


def propsort(data,pn):
    '''sort datasets in data using the property named pn as key'''
    data.sort(cmp=lambda x,y:cmp(eval(x[0].props[pn]),eval(y[0].props[pn])))
    
listobs=['Green_'+str(flavor)]
ll=pyalps.load.Hdf5Loader()
result_files = pyalps.getResultFiles(prefix=res_file)

data_G = ll.ReadDMFTIterations(result_files, 'G_omega', measurements=listobs, verbose=True)
data_G0 = ll.ReadDMFTIterations(result_files, 'G0_omega', measurements=listobs, verbose=True)
grouped_G = pyalps.groupSets(pyalps.flatten(data_G), ['iteration'])   # [iteration][result_files(assumed only 1) x measurements(here only 1)]
grouped_G0 = pyalps.groupSets(pyalps.flatten(data_G0), ['iteration'])   # [iteration][result_files(assumed only 1) x measurements(here only 1)]
propsort(grouped_G,'iteration')  # however, the order may be messed up in the function pyalps.plot.plot()
propsort(grouped_G0,'iteration')
nd=[]
for i in range(0,len(grouped_G)):  # assumed: len(grouped_G)==len(grouped_G0)
    G = grouped_G[i][0]
    G0 = grouped_G0[i][0]
    r = pyalps.DataSet()
    r.y = np.zeros(len(G.y))
    for w in range(0,len(G.y)):
      selfenergy = 1./G0.y[w] - 1./G.y[w]
      if real_part:
        r.y[w] = selfenergy.real
      else:
        r.y[w] = selfenergy.imag
      
    r.x = np.array([(2.*e+1)*pi/G.props['BETA'] for e in G.x])
    r.props = G.props
    r.props['label'] = "it"+r.props['iteration']
    r.props['line']="scatter"
    r.props['fillmarkers'] = False
    nd.append( r )

plt.figure()
plt.xlabel(r'$i\omega_n$')
if real_part:
  plt.ylabel(r'$Re \Sigma_{flavor='+str(flavor)+r'}(i\omega_n)$')
else:
  plt.ylabel(r'$Im \Sigma_{flavor='+str(flavor)+r'}(i\omega_n)$')  
plt.title(result_files[0])
pyalps.plot.plot(nd)
plt.legend()
plt.show()
