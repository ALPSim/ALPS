# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2015 by Synge Todo <wistaria@comp-phys.org> 
# 
# SPDX-License-Identifier: MIT
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
