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

basename = 'parm_complex'

#prepare the input parameters
parms = []
parms.append( { 
    'LATTICE'                               : "open ladder",
    'L'                                     : 10,
    'MODEL_LIBRARY'                         : "mymodels.xml",
    'MODEL'                                 : "fermion Hubbard",
    'CONSERVED_QUANTUMNUMBERS'              : 'Nup,Ndown',
    'Nup_total'                             : 10,
    'Ndown_total'                           : 10,
    't0'                                    : "1+0.6*I",
    'ct0'                                   : "1-0.6*I",
    't1'                                    : 0.1,
    'U'                                     : 0.,
    'SWEEPS'                                : 6,
    'MAXSTATES'                             : 400,
    'COMPLEX'                               : 1,
   } )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles(basename,parms)
res = pyalps.runApplication('mps_optim',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix=basename), ['Energy'])

en_exact = -28.1129977
print 'Exact energy for MAXSTATES=inf ::', en_exact
for d in pyalps.flatten(data):
    print d.props['observable'], '=', d.y

