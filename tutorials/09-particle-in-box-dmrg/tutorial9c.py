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

#prepare the input parameters
parms = { 
          'LATTICE'          : "next-nearest chain lattice", 
          'LATTICE_LIBRARY'   : "lattices_dmrg.xml", 
          'MODEL'             : "spin",
          'L'                 : 20,
          't'                 : 1,
          't1'                : 0.01,
          'K'                 : '2*3.1415927*2/L',
          'V'                 : 'cos(K*x)+cos(3*K*x)/9+cos(5*K*x)/25+cos(7*K*x)/36-2',
          'SWEEPS'            : 100,
          'WAVEFUNCTION_FILE' : "psi2.dat",
          'OUTPUT_LEVEL'      : 1
        }

#write the input file and run the simulation

input_file = pyalps.writeParameterFile('parm9c',parms)
res = pyalps.runApplication('simple_dmrg',input_file)

# load the text file and plot the wave function
raw = np.loadtxt(parms['WAVEFUNCTION_FILE']).transpose()
data = pyalps.DataSet()
data.x = raw[0]
data.y = raw[1]
data.props['xlabel'] = 'x'

plt.figure()
pyalps.pyplot.plot(data)
plt.show()
