# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2025 by ALPS Collaboration
#
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.plot
import pyalps.fit_wrapper as fw

#prepare the input parameters
parms= []
for lattice in [32, 64, 96, 128]:
    parms.append({
            'LATTICE'                   : "open chain lattice",
            'MODEL'                     : "spin",
            'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
            'Sz_total'                  : 0,
            'J'                         : 1,
            'SWEEPS'                    : 4,
            'L'                         : lattice,
            'MAXSTATES'                 : 100,
            'NUMBER_EIGENVALUES'        : 2
        })

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm_spin_one_half_gap_multiple',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load all measurements for all lattices
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm_spin_one_half_gap_multiple'))

#sort results by the lattice length
sorted_data = sorted(data, key=lambda x: x[0].props['L'])

# create data set for plot: gap vs. 1/L
gapplot = pyalps.DataSet()
gapplot.props['xlabel']='$1/L$'
gapplot.props['ylabel']='Gap $\Delta/J$'
gapplot.props['label']='D=100'
gapplot.props['line']='.'

x = []
y = []
for measure in sorted_data:
    for s in measure:
        if s.props['observable'] == 'Energy':
            L = s.props['L']
            iL = 1.0/L
            gap = abs(s.y[1] - s.y[0])
            s.props['gap'] = gap
            x.append(iL)
            y.append(gap)

gapplot.x = x
gapplot.y = y

# plot the gap vs. 1/L curve:
plt.figure()
pyalps.plot.plot(gapplot)
plt.legend()
plt.xlim(0,0.04)
plt.ylim(0,0.2)

# fit the curve with a linear function
pars = [fw.Parameter(0.1), fw.Parameter(0.2)]
f = lambda self, x, p: p[0]()+p[1]()*x
fw.fit(None, f, pars, np.array(gapplot.y), np.array(gapplot.x))

# plot the fitted curve
x = np.linspace(0.0, 0.035, 100)
plt.plot(x, f(None,x,pars))

print("Gap at thermodynamic limit: ", pars[0]())

plt.show()

