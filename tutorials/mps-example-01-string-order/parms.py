from __future__ import print_function
# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2016 by Michele Dolfi <dolfim@phys.ethz.ch>
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

basename = 'parm_string'

#prepare the input parameters
parms = []
parms.append( { 
    'LATTICE'                               : "inhomogeneous chain lattice",
    'L'                                     : 12,
    'MODEL_LIBRARY'                         : "mymodels.xml",
    'MODEL'                                 : "boson Hubbard",
    'Nmax'                                  : 4,
    'CONSERVED_QUANTUMNUMBERS'              : 'N',
    'N_total'                               : 6,
    't'                                     : 1,
    'U'                                     : 8.,
    'SWEEPS'                                : 5,
    'NUMBER_EIGENVALUES'                    : 1,
    'MAXSTATES'                             : 100,
    'MEASURE_LOCAL[Local density]'          : 'n',
    'MEASURE_LOCAL_AT[String order 2]'      : 'st:st|(4,5),(5,6),(6,7)',
    'MEASURE_LOCAL_AT[String order 4]'      : 'st:st:st:st|((4,5,6,7),(3,4,5,6),(5,6,7,8))',
   } )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles(basename,parms)
res = pyalps.runApplication('mps_optim',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix=basename), ['String order 2', 'String order 4'])

for d in pyalps.flatten(data):
    print('##', d.props['observable'])
    for x,y in zip(d.x, d.y[0]):
        print('Sites:', x)
        print('Val:  ', y)

