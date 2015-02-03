# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Jan Gukelberger <gukelberger@phys.ethz.ch> 
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

#prepare the input parameters
parms=[]
for l in [6, 8, 10]:
    parms.append(
      { 
        'LATTICE'                   : "ladder", 
        'MODEL'                     : "spin",
        'local_S'                   : 0.5,
        'J0'                        : 1,
        'J1'                        : 1,
        'L'                         : l,
        'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
        'Sz_total'                  : 0
      }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('ed03b',parms)
res = pyalps.runApplication('sparsediag',input_file)

#load all measurements for all states
data = pyalps.loadSpectra(pyalps.getResultFiles(prefix='ed03b'))

# collect spectra over all momenta for every simulation
spectra = {}
for sim in data:
  l = int(sim[0].props['L'])
  all_energies = []
  spectrum = pyalps.DataSet()
  for sec in sim:
    all_energies += list(sec.y)
    spectrum.x = np.concatenate((spectrum.x,np.array([sec.props['TOTAL_MOMENTUM'] for i in range(len(sec.y))])))
    spectrum.y = np.concatenate((spectrum.y,sec.y))
  spectrum.y -= np.min(all_energies)
  spectrum.props['line'] = 'scatter'
  spectrum.props['label'] = 'L='+str(l)
  spectra[l] = spectrum


# plot
plt.figure()
pyalps.plot.plot(spectra.values())
plt.legend()
plt.title('Antiferromagnetic Heisenberg ladder (S=1/2)')
plt.ylabel('Energy')
plt.xlabel('Momentum')
plt.xlim(0,2*3.1416)
plt.ylim(0,2.5)
plt.show()
