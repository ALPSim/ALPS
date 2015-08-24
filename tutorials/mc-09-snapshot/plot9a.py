# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2015 by Synge Todo <wistaria@comp-phys.org> 
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
import pyalps.plot as alpsplot
import matplotlib.pyplot as pyplot

data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm9a'),
    ['Specific Heat', 'Magnetization Density^2', 'Binder Ratio of Magnetization'])
for item in pyalps.flatten(data):
    item.props['L'] = int(item.props['L'])

magnetization2 = pyalps.collectXY(data, x='T', y='Magnetization Density^2', foreach=['L'])
magnetization2.sort(key=lambda item: item.props['L'])

specificheat = pyalps.collectXY(data, x='T', y='Specific Heat', foreach=['L'])
specificheat.sort(key=lambda item: item.props['L'])

binderratio = pyalps.collectXY(data, x='T', y='Binder Ratio of Magnetization', foreach=['L'])
binderratio.sort(key=lambda item: item.props['L'])

pyplot.figure()
alpsplot.plot(magnetization2)
pyplot.xlabel('Temperture $T$')
pyplot.ylabel('Magnetization Density Squared $m^2$')
pyplot.legend(loc='best')

pyplot.figure()
alpsplot.plot(specificheat)
pyplot.xlabel('Temperature $T$')
pyplot.ylabel('Specific Heat C')
pyplot.legend(loc='best')

pyplot.figure()
alpsplot.plot(binderratio)
pyplot.xlabel('Temperature $T$')
pyplot.ylabel('Binder Ratio of Magnetization')
pyplot.legend(loc='best')

pyplot.show()
