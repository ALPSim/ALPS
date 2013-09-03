#############################################################################
#
# ALPS Project Applications: Directed Worm Algorithm  
#
# Copyright (C) 2013 by Lode Pollet      <pollet@phys.ethz.ch>  
#                       Ping Nang Ma     <pingnang@phys.ethz.ch> 
#                       Matthias Troyer  <troyer@phys.ethz.ch>    
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
#############################################################################

# The headers
import pyalps

# Set up a python list of parameters (python) dictionaries:
parms = [{
  'LATTICE'         : "square lattice",          
  'MODEL'           : "boson Hubbard",
  'L'               : 20,
  'Nmax'            : 20,
  't'               : 1.,
  'U'               : 16.,
  'mu'              : 32.,
  'T'               : 1.,
  'THERMALIZATION'  : 10000,
  'SWEEPS'          : 100000,
  'SKIP'            : 400
}]

# Write into XML input file:
input_file = pyalps.writeInputFiles('parm1a',parms)

# and run the application dwa:
pyalps.runApplication('dwa', input_file, Tmin=10, writexml=True)

# We first get the list of all hdf5 result files via:
files = pyalps.getResultFiles(prefix='parm1a', format='hdf5')

# and then extract, say the timeseries of the Density measurements:
ar = pyalps.hdf5.h5ar(files[0])
density_timeseries = ar['/simulation/results']['Density']['timeseries']['data']

# We can then visualize graphically:
import matplotlib.pyplot as plt
plt.plot(density_timeseries)
plt.show()

# ALPS Python provides a convenient tool to check whether a measurement observable(s) has (have) reached steady state equilibrium.
#
# Here is one example:
pyalps.steady_state_check(files[0], 'Density')

# and another one:
pyalps.steady_state_check(files[0], ['Density', 'Energy Density'])

# To see the complete log
pyalps.steady_state_check(files[0], ['Density', 'Energy Density'], includeLog=True)
