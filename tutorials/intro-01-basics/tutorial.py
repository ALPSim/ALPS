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
import pyalps.pyplot

#prepare the input parameters
parms = []
for t in [1.5,2,2.5]:
    parms.append(
        { 
          'LATTICE'        : "square lattice", 
          'T'              : t,
          'J'              : 1 ,
          'THERMALIZATION' : 1000,
          'SWEEPS'         : 100000,
          'UPDATE'         : "cluster",
          'MODEL'          : "Ising",
          'L'              : 8
        }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm1',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5,writexml=True)

#get the list of result files
result_files = pyalps.getResultFiles(prefix='parm1')
print result_files

#print the observables stored in those files:
print pyalps.loadObservableList(result_files)

#load a selection of measurements:
data = pyalps.loadMeasurements(result_files,['|Magnetization|','Magnetization^2'])

#calculate the Binder cumulants




#make one plot with all data
for dataset in binning:
    dataset.props['label'] = 'L='+str(dataset.props['L'])

plt.figure()
plt.title('Binning analysis for cluster updates')
plt.xlabel('binning level')
plt.ylabel('Error of |Magnetization|')
pyalps.pyplot.plot(binning)
plt.legend()
plt.show()


# make individual plots for each system size
for dataset in binning:
    plt.figure()
    plt.title('Binning analysis for L='+str(dataset.props['L']))
    plt.xlabel('binning level')
    plt.ylabel('Error of |Magnetization|')
    pyalps.pyplot.plot(dataset)

plt.show()
