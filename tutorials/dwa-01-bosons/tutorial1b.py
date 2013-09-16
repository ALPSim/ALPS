# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2013 by Matthias Troyer <troyer@phys.ethz.ch>,
#                       Ping Nang Ma    <pingnang@phys.ethz.ch> 
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


# Preparing and running the simulation using Python
import pyalps

parms = []
for L in [4,6,8]:
  for t in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
    parms.append(
        { 
          'LATTICE'                 : "square lattice", 
          'MODEL'                   : "boson Hubbard",
          'T'                       : 0.1,
          'L'                       : L ,
          't'                       : t ,
          'mu'                      : 0.5,
          'U'                       : 1.0 ,
          'Nmax'                    : 2 ,
          'THERMALIZATION'          : 100000,
          'SWEEPS'                  : 2000000,
          'SKIP'                    : 500,
          'MEASURE[Winding Number]': 1
        }
    )

input_file = pyalps.writeInputFiles('parm1b', parms)
res = pyalps.runApplication('dwa', input_file, Tmin=5, writexml=True)


# Evaluating the simulation and preparing plots using Python
import pyalps
import matplotlib.pyplot as plt
import pyalps.plot as aplt

data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm1b'),'Stiffness')
rhos = pyalps.collectXY(data,x='t',y='Stiffness',foreach=['L'])

for rho in rhos:
  rho.y = rho.y * float(rho.props['L'])

plt.figure()
aplt.plot(rhos)
plt.xlabel('Hopping $t/U$')
plt.ylabel('$\\rho _sL$')
plt.legend()
plt.title('Scaling plot for Bose-Hubbard model')
plt.show()

