# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Matthias Troyer <troyer@phys.ethz.ch> 
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

import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.plot

data = pyalps.loadSpectra(pyalps.getResultFiles(prefix='parm2a'))
data += pyalps.loadSpectra(pyalps.getResultFiles(prefix='parm2b'))

lengths = []
min_energies = {}

# extract the ground state energies over all momenta for every simulation
for sim in data:
  l = int(sim[0].props['L'])
  if l not in lengths: lengths.append(l)
  sz = int(sim[0].props['Sz_total'])
  s = float(sim[0].props['local_S'])
  all_energies = []
  for sec in sim:
    all_energies += list(sec.y)
  min_energies[(l,s,sz)]= np.min(all_energies)
  

# make a plot of the triplet gap as function of system size 
plt.figure()
for s in [0.5,1]:
  gapplot = pyalps.DataSet()
  gapplot.x = 1./np.sort(lengths)
  gapplot.y = [min_energies[(l,s,1)] -min_energies[(l,s,0)] for l in np.sort(lengths)]  
  gapplot.props['xlabel']='$1/L$'
  gapplot.props['ylabel']='Triplet gap $\Delta/J$'
  gapplot.props['label']='S='+str(s)
  pyalps.plot.plot(gapplot)

plt.legend()
plt.xlim(0,0.25)
plt.ylim(0,1.0)
plt.show()
