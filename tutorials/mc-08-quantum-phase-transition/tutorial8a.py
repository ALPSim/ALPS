from __future__ import print_function
# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Brigitte Surer
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
import pyalps.plot
import numpy as np
import pyalps.fit_wrapper as fw
from math import sqrt

#prepare the input parameters
parms = []
for j2 in [0.,1.]:
    for t in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0]:
        parms.append(
            { 
              'LATTICE'        : "coupled ladders", 
              'local_S'        : 0.5,
              'ALGORITHM'      : 'loop',
              'SEED'           : 0,
              'T'              : t,
              'J0'             : 1 ,
              'J1'             : 1,
              'J2'             : j2,
              'THERMALIZATION' : 5000,
              'SWEEPS'         : 50000, 
              'MODEL'          : "spin",
              'L'              : 8,
              'W'              : 4
            }
    )
    
#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm8a',parms)
pyalps.runApplication('loop',input_file)

data = pyalps.loadMeasurements(pyalps.getResultFiles(pattern='parm8a.task*.out.h5'),['Staggered Susceptibility','Susceptibility'])
susc1=pyalps.collectXY(data,x='T',y='Susceptibility', foreach=['J2'])

lines = []
for data in susc1:
    pars = [fw.Parameter(1), fw.Parameter(1)]
    data.y= data.y[data.x < 1]
    data.x= data.x[data.x < 1]
    f = lambda self, x, pars: (pars[0]()/np.sqrt(x))*np.exp(-pars[1]()/x)
    fw.fit(None, f, pars, np.array([v.mean for v in data.y]), data.x)
    prefactor = pars[0].get()
    gap = pars[1].get()
    print(prefactor,gap)
    
    lines += plt.plot(data.x, f(None, data.x, pars))
    lines[-1].set_label('$J_2=%.4s$: $\chi = \\frac{%.4s}{T}\exp(\\frac{-%.4s}{T})$' % (data.props['J2'], prefactor,gap))

plt.figure()
pyalps.plot.plot(susc1)
plt.xlabel(r'$T$')
plt.ylabel(r'$\chi$')
plt.title('gap is %.4s' % gap)
plt.legend()
plt.show()
