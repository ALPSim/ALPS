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
import pyalps.pyplot

#prepare the input parameters

parms=[]
for sz in [0, 1]:
      parms.append(
        {
          'LATTICE'                   : "ladder", 
          'MODEL'                     : "spin",
          'local_S'                   : 0.5,
          'J0'                        : 1,
          'J1'                        : 1,
          'L'                         : 8,
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'Sz_total'                  : sz,
          'NUMBER_EIGENVALUES'        : 10
        }
      )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm1e',parms)
res = pyalps.runApplication('sparsediag',input_file)

#load all measurements for all states
data = pyalps.loadSpectra(pyalps.getResultFiles(prefix='parm1e'))

energies=[]
# get ground state energy
for s in pyalps.flatten(data):
    energies += list(s.y)

groundstate_energy = np.min(energies)

#plot spectra
spectrumplot = []
for set in data:
  plot = pyalps.DataSet()
  plot.props['label']='Sz='+str(set[0].props['Sz_total'])
  plot.props['line'] = 'scatter'
  for s in set:
    plot.x = np.concatenate((plot.x,np.array([s.props['TOTAL_MOMENTUM'] for i in range(0,len(s.y))])))
    plot.y = np.concatenate((plot.y,s.y - groundstate_energy))
  spectrumplot.append(plot)
  
spectrumplot[0] = pyalps.subtract_spectrum(spectrumplot[0],spectrumplot[1],tolerance=1e-12)

plt.figure()
pyalps.pyplot.plot(spectrumplot)
plt.legend()
plt.xlim(0,2*3.1416)
plt.ylabel('Energy')
plt.ylim(0,3)
plt.show()
