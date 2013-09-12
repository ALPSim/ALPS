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
# WARNING: With python library scripts of revision 6194 and older you get all markers in black.
#



import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.plot
from math import pi

print "DESCRIPTION: This script is designed for comparison of final G, G0 and selfenergy for all result files in the current directory (recursively). The selfenergy is calculated via Dyson equation from the stored G and G0."
print

result_files_ = pyalps.getResultFiles()
result_files = []
print "RESULT FILES to be processed :"
for a in result_files_:
    res_file = [a]
    obs = pyalps.loadObservableList(res_file)
    if obs[0].count('G_tau')>0:
      print "  ",a
      result_files.append( a )

print "Please select flavor [default: 0] :"
answer = raw_input('--> ')
flavor = 0
if len(answer)>0:
  flavor = eval(answer)
listobs=[str(flavor)]   # previous format: "Green_"+str(flavor)

data_G_tau = pyalps.loadMeasurements(result_files, respath='/simulation/results/G_tau', measurements=listobs, verbose=True)  # [result_files][measurements(here only 1)]
for d in pyalps.flatten(data_G_tau):
    d.x = d.x*d.props["BETA"]/float(d.props["N"])
    d.props['line']="scatter"
    d.props['fillmarkers'] = False
plt.figure()
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G_{flavor='+str(flavor)+r'}(\tau)$')
pyalps.plot.plot(data_G_tau)
plt.legend()

data_G0_tau = pyalps.loadMeasurements(result_files, respath='/simulation/results/G0_tau', measurements=listobs, verbose=True)
for d in pyalps.flatten(data_G0_tau):
    d.x = d.x*d.props["BETA"]/float(d.props["N"])
    d.props['line']="scatter"
    d.props['fillmarkers'] = False
plt.figure()
plt.xlabel(r'$\tau$')
plt.ylabel(r'$G^0_{flavor='+str(flavor)+r'}(\tau)$')
pyalps.plot.plot(data_G0_tau)
plt.legend()

data_G_omega = pyalps.loadMeasurements(result_files, respath='/simulation/results/G_omega', measurements=listobs, verbose=True)
data_G_omegareal = []
for d in pyalps.flatten(data_G_omega):
    r = pyalps.DataSet()
    r.y = np.array(d.y.real)
    r.x = np.array([(2.*e+1)*pi/d.props['BETA'] for e in d.x])
    r.props = d.props
    r.props['line']="scatter"
    r.props['fillmarkers'] = False
    data_G_omegareal.append( r )
plt.figure()
plt.xlabel(r'$i\omega_n$')
plt.ylabel(r'$Re G_{flavor='+str(flavor)+r'}(i\omega_n)$')
pyalps.plot.plot(data_G_omegareal)
plt.legend()
data_G_omegaimag = []
for d in pyalps.flatten(data_G_omega):
    r = pyalps.DataSet()
    r.y = np.array(d.y.imag)
    r.x = np.array([(2.*e+1)*pi/d.props['BETA'] for e in d.x])
    r.props = d.props
    r.props['line']="scatter"
    r.props['fillmarkers'] = False
    data_G_omegaimag.append( r )
plt.figure()
plt.xlabel(r'$i\omega_n$')
plt.ylabel(r'$Im G_{flavor='+str(flavor)+r'}(i\omega_n)$')
pyalps.plot.plot(data_G_omegaimag)
plt.legend()

data_G0_omega = pyalps.loadMeasurements(result_files, respath='/simulation/results/G0_omega', measurements=listobs, verbose=True)
data_G0_omegareal = []
for d in pyalps.flatten(data_G0_omega):
    r = pyalps.DataSet()
    r.y = np.array(d.y.real)
    r.x = np.array([(2.*e+1)*pi/d.props['BETA'] for e in d.x])
    r.props = d.props
    r.props['line']="scatter"
    r.props['fillmarkers'] = False
    data_G0_omegareal.append( r )
plt.figure()
plt.xlabel(r'$i\omega_n$')
plt.ylabel(r'$Re G^0_{flavor='+str(flavor)+r'}(i\omega_n)$')
pyalps.plot.plot(data_G0_omegareal)
plt.legend()
data_G0_omegaimag = []
for d in pyalps.flatten(data_G0_omega):
    r = pyalps.DataSet()
    r.y = np.array(d.y.imag)
    r.x = np.array([(2.*e+1)*pi/d.props['BETA'] for e in d.x])
    r.props = d.props
    r.props['line']="scatter"
    r.props['fillmarkers'] = False
    data_G0_omegaimag.append( r )
plt.figure()
plt.xlabel(r'$i\omega_n$')
plt.ylabel(r'$Im G^0_{flavor='+str(flavor)+r'}(i\omega_n)$')
pyalps.plot.plot(data_G0_omegaimag)
plt.legend()

selfenergy_real=[]
for i in range(0,len(data_G0_omega)):
    G = data_G_omega[i][0]
    G0 = data_G0_omega[i][0]
    r = pyalps.DataSet()
    r.y = np.zeros(len(G.y))
    for w in range(0,len(G.y)):
      selfenergy = 1./G0.y[w] - 1./G.y[w]
      r.y[w] = selfenergy.real    
    r.x = np.array([(2.*e+1)*pi/G.props['BETA'] for e in G.x])
    r.props = G.props
    r.props['line']="scatter"
    r.props['fillmarkers'] = False
    selfenergy_real.append( r )
plt.figure()
plt.xlabel(r'$i\omega_n$')
plt.ylabel(r'$Re \Sigma_{flavor='+str(flavor)+r'}(i\omega_n)$')
pyalps.plot.plot(selfenergy_real)
plt.legend()
selfenergy_imag=[]
for i in range(0,len(data_G0_omega)):
    G = data_G_omega[i][0]
    G0 = data_G0_omega[i][0]
    r = pyalps.DataSet()
    r.y = np.zeros(len(G.y))
    for w in range(0,len(G.y)):
      selfenergy = 1./G0.y[w] - 1./G.y[w]
      r.y[w] = selfenergy.imag
    r.x = np.array([(2.*e+1)*pi/G.props['BETA'] for e in G.x])
    r.props = G.props
    r.props['line']="scatter"
    r.props['fillmarkers'] = False
    selfenergy_imag.append( r )
plt.figure()
plt.xlabel(r'$i\omega_n$')
plt.ylabel(r'$Im \Sigma_{flavor='+str(flavor)+r'}(i\omega_n)$')
pyalps.plot.plot(selfenergy_imag)
plt.legend()

plt.show()

