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
for t in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
    parms.append(
        { 
          'LATTICE'        : "square lattice", 
          'MODEL'          : "boson Hubbard",
          'T'              : 0.1,
          'L'              : 4 ,
          't'              : t ,
          'mu'             : 0.5,
          'U'              : 1.0 ,
          'NONLOCAL'       : 0 ,
          'Nmax'           : 2 ,
          'THERMALIZATION' : 10000,
          'SWEEPS'         : 500000
        }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('mc05a',parms)
res = pyalps.runApplication('worm',input_file,Tmin=5)

#load the magnetization and collect it as function of field h
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='mc05a'),'Stiffness')
rhos = pyalps.collectXY(data,x='t',y='Stiffness')

#make plot
plt.figure()
pyalps.plot.plot(rhos)
plt.xlabel('Hopping $t/U$')
plt.ylabel('Superfluid density $\ho _s$')
plt.show()
