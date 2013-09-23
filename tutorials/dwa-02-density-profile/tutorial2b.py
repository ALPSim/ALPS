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

parms = [
  {
    'LATTICE' : 'inhomogeneous simple cubic lattice' ,
    'L'       : 60 ,

    'MODEL'   : 'boson Hubbard' ,
    'Nmax'    : 20 ,

    't'  : 1. ,
    'U'  : 60. ,
    'mu' : '40. - (0.09416*(x-(L-1)/2.)*(x-(L-1)/2.) + 0.12955*(y-(L-1)/2.)*(y-(L-1)/2.) + 0.11496*(z-(L-1)/2.)*(z-(L-1)/2.))' ,

    'T'  : 1. ,

    'THERMALIZATION' : 1000000 ,
    'SWEEPS'         : 3000000 ,
    'SKIP'           : 1000 , 

    'MEASURE[Local Density]': 1
  }
]

input_file = pyalps.writeInputFiles('parm2b', parms)
res = pyalps.runApplication('dwa', input_file)


# Evaluating and plotting in Python
import pyalps
import pyalps.plot as aplt;

data     = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm2b'), 'Local Density');
aplt.plot3D(data, centeredAtOrigin=True, layer="center")

