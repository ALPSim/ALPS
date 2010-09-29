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

import pyalps
import matplotlib.pyplot as plt
import pyalps.plot

#prepare the input parameters
parms = []
for t in [0.045, 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.065]:
  for l in [4, 6, 8]:
    parms.append(
        { 
          'LATTICE'        : "square lattice", 
          'MODEL'          : "boson Hubbard",
          'T'              : 0.05,
          'L'              : l ,
          't'              : t ,
          'U'              : 1.0 ,
          'mu'             : 0.5, 
          'NONLOCAL'       : 0 ,
          'Nmax'           : 2 ,
          'THERMALIZATION' : 15000,
          'SWEEPS'         : 600000
        }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm5b',parms)
res = pyalps.runApplication('worm',input_file,Tmin=5)

#load the magnetization and collect it as function of field h
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm5b'),'Stiffness')
rhos = pyalps.collectXY(data,x='t',y='Stiffness',foreach=['L'])

# multiply with the system size for the scaling plot
for s in rhos:
  s.y = s.y * float(s.props['L'])

#make plot
plt.figure()
pyalps.plot.plot(rhos)
plt.xlabel('Hopping $t/U$')
plt.ylabel('$\\rho _sL$')
plt.legend()
plt.title('Scaling plot for Bose-Hubbard model')
plt.show()
