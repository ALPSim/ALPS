# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Matthias Troyer <troyer@phys.ethz.ch> 
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

# Please run all four other tutorials before running this one. 
# This tutorial relies on the results created in those tutorials

import pyalps
import matplotlib.pyplot as plt
import pyalps.plot

# load all files
data = pyalps.loadMeasurements(pyalps.getResultFiles(),'Susceptibility')

#flatten the hierarchical structure
data = pyalps.flatten(data)

# collect the susceptibility
susceptibility = pyalps.collectXY(data,x='T',y='Susceptibility',foreach=['MODEL','LATTICE'])

# assign labels to the data depending on the properties
for s in susceptibility:
  # print s.props
  if s.props['LATTICE']=='chain lattice':
    s.props['label'] = "chain"
  elif s.props['LATTICE']=='ladder':
    s.props['label'] = "ladder"
  if s.props['MODEL']=='spin':
    s.props['label'] = "quantum " + s.props['label']
  elif s.props['MODEL']=='Heisenberg':
    s.props['label'] = "classical " + s.props['label']

#make plot
plt.figure()
pyalps.plot.plot(susceptibility)
plt.xlabel('Temperature $T/J$')
plt.ylabel('Susceptibility $\chi J$')
plt.ylim(0,0.25)
plt.legend()
plt.show()

