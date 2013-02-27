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

print "DESCRIPTION: This script gets the iteration-resolved observable Green's function from the specifyied result file. For single site problems only."
print

all_result_files = pyalps.getResultFiles()
if len(all_result_files)>1:
  print "Choose the (single) result file, which is to be examined: "
  for i in range(0,len(all_result_files)):
    print i,':', all_result_files[i]
  result_files = [all_result_files[eval(raw_input('--> '))]]
else:
  print "Result file to be examined: "
  print all_result_files[0]
  result_files=all_result_files

#print "Enter the (prefix of the) result file(s), which is(are) to examine :"
#res_file = raw_input('--> ')
#result_files = pyalps.getResultFiles(prefix=res_file)

print "Do you want to get the Green's function in the Matsubara frequency [default] or in the imaginary-time [type: 't'] representation ?"
answer1 = raw_input('--> ')
tau_repr = False
if len(answer1)>0 and answer1[0]=='t':
  tau_repr = True

real_part = True
if not tau_repr:
  print "Do you want to plot the imaginary part of the Green's function [default] or the real part [type: 'r'] ?"
  answer2 = raw_input('--> ')
  if len(answer2)==0 or answer2[0]!='r':
    real_part = False

print "Do you want to get the full impurity Green's function [default] or the Weiss field (G^0) [type: '0'] ?"
answer3 = raw_input('--> ')
G0 = False
if len(answer3)>0 and answer3[0]=='0':
  G0 = True

print "Please select flavor [default: 0] :"
answer4 = raw_input('--> ')
flavor = 0
if len(answer4)>0:
  flavor = eval(answer4)


obs='G'
if G0:
  obs += '0'
if tau_repr:
  obs += '_tau'
else:
  obs += '_omega'

def propsort(data,pn):
    '''sort datasets in data using the property named pn as key'''
    data.sort(cmp=lambda x,y:cmp(eval(x[0].props[pn]),eval(y[0].props[pn])))
    
listobs=[str(flavor)]  # previous format: "Green_"+str(flavor)
ll=pyalps.load.Hdf5Loader()
data = ll.ReadDMFTIterations(result_files, observable=obs, measurements=listobs, verbose=True)
grouped = pyalps.groupSets(pyalps.flatten(data), ['iteration'])   # [iteration][result_files x measurements(here only 1)]
propsort(grouped,'iteration')  # however, the order is messed up in the function pyalps.plot.plot()
nd=[]
for i in range(0,len(grouped)):
  for dat_set in grouped[i]:
    r = pyalps.DataSet()
    if real_part:
      r.y = np.array(dat_set.y.real)
    else:
      r.y = np.array(dat_set.y.imag)
    if tau_repr:
      r.x = np.array([e*dat_set.props['BETA']/float(dat_set.props['N']) for e in dat_set.x])
    else:
      r.x = np.array([(2.*e+1)*pi/dat_set.props['BETA'] for e in dat_set.x])
    r.props = dat_set.props
    r.props['label'] = "it"+r.props['iteration']
    r.props['line']="scatter"
    r.props['fillmarkers'] = False
    nd.append( r )

plt.figure()
xlab=''
if tau_repr:
  xlab=r'$\tau$'
else:
  xlab=r'$i\omega_n$'
ylab=''
if G0:
  ylab=r'$G^0$'
else:
  if real_part:
    ylab=r'$Re G$'
  else:
    ylab=r'$Im G$'

plt.xlabel(xlab)
plt.ylabel(ylab+r'$_{flavor='+str(flavor)+r'}($'+xlab+r'$)$')
title_=''
for i in range(0,len(result_files)):
  if i>0:
    title_+=", "
  title_+=result_files[i]
plt.title(title_)
pyalps.plot.plot(nd)
plt.legend()
plt.show()
