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
import matplotlib.pyplot as plt
import pyalps.pyplot
import numpy as np

#prepare the input parameters
parms = [{ 
          'LATTICE'                   : "ladder", 
          'MODEL'                     : "spin",
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'local_S'                   : 0.5,
          'J0'                        : 1,
          'J1'                        : 1,
          'L'                         : 6
        }]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm2b',parms)
res = pyalps.runApplication('fulldiag',input_file)

#run the evaluation and load all the plots
data = pyalps.evaluateFulldiagVersusT(pyalps.getResultFiles(prefix='parm2b'),DELTA_T=0.05, T_MIN=0.05, T_MAX=5.0)

#make plot
for s in pyalps.flatten(data):
  plt.figure()
  plt.title("Antiferromagnetic Heisenberg chain")
  pyalps.pyplot.plot(s)


# make a plot of the spectrum, first load all measurements for all states
spectrum = pyalps.loadSpectra(pyalps.getResultFiles(prefix='parm2b'))

energies=[]
# get ground state energy
for set in spectrum:
  for s in set:
    energies += list(s.y)

groundstate_energy = np.min(energies)

#collect spctra for each Sz into one DataSet in a Python dictionary
spectrumplot = {}
for sz in [0,1,2]:
    spectrumplot[sz] = pyalps.DataSet()
    spectrumplot[sz].props['label']='Sz='+str(sz)
    spectrumplot[sz].props['line'] = 'scatter'

for s in spectrum[0]:
  sz = int(s.props['Sz'])
  if sz in spectrumplot:
    spectrumplot[sz].x = np.concatenate((spectrumplot[sz].x,np.array([s.props['TOTAL_MOMENTUM'] for i in range(0,len(s.y))])))
    spectrumplot[sz].y = np.concatenate((spectrumplot[sz].y,s.y - groundstate_energy))

plt.figure()
pyalps.pyplot.plot(spectrumplot.values())
plt.legend()
plt.xlim(0,2*3.1416)
plt.ylabel('Energy')
plt.ylim(0,3)
plt.title('Antiferromagnetic Heisenberg ladder')
plt.show()


