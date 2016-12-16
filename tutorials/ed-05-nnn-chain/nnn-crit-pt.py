from __future__ import print_function
# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Bela Bauer <bauerb@phys.ethz.ch> 
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
import pyalps.plot
import numpy as np
import matplotlib.pyplot as plt
import copy
import math
from pyalps import flatten

prefix = 'alps-nnn-heisenberg'
parms = []
for L in [6,8]:
    for Szt in [0,1]:
      for J1 in np.linspace(0,0.5,6):
          parms.append({
              'LATTICE'              : "nnn chain lattice",
              'MODEL'                : "spin",
              'local_S'              : 0.5,
              'J'                    : 1,
              'NUMBER_EIGENVALUES'   : 2,
              'CONSERVED_QUANTUMNUMBER' : 'Sz',
              'Sz_total'             : Szt,
              'J1'                   : J1,
              'L'                    : L
          })

input_file = pyalps.writeInputFiles(prefix,parms)
res = pyalps.runApplication('sparsediag', input_file)
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix=prefix))

# join all momenta
grouped = pyalps.groupSets(pyalps.flatten(data), ['J1', 'L', 'Sz_total'])
nd = []
for group in grouped:
    ally = []
    allx = []
    for q in group:
        ally += list(q.y)
        allx += list(q.x)
    r = pyalps.DataSet()
    sel = np.argsort(ally)
    r.y = np.array(ally)[sel]
    r.x = np.array(allx)[sel]
    r.props = pyalps.dict_intersect([q.props for q in group])
    nd.append( r )
data = nd

# remove states in the s=1 sector from the s=0 sector
grouped = pyalps.groupSets(pyalps.flatten(data), ['J1', 'L'])
nd = []
for group in grouped:
    if group[0].props['Sz_total'] == 0:
        s0 = group[0]
        s1 = group[1]
    else:
        s0 = group[1]
        s1 = group[0]
    s0 = pyalps.subtract_spectrum(s0, s1)
    nd.append(s0)
    nd.append(s1)
data = nd

# make a plot of the ground state ('gs') and the first excited state ('fe')
# as a function of J_1
sector_E = []
grouped = pyalps.groupSets(pyalps.flatten(data), ['Sz_total', 'J1', 'L'])
for group in grouped:
    allE = []
    for q in group:
        allE += list(q.y)
    allE = np.sort(allE)
    
    d = pyalps.DataSet()
    d.props = pyalps.dict_intersect([q.props for q in group])
    d.x = np.array([0])
    d.y = np.array([allE[0]])
    d.props['which'] = 'gs'
    d.props['line'] = '.-'
    sector_E.append(d)
    
    d2 = copy.deepcopy(d)
    d2.y = np.array([allE[1]])
    d2.props['which'] = 'fe'
    d2.props['line'] = '.-'
    sector_E.append(d2)

sector_energies = pyalps.collectXY(sector_E, 'J1', 'Energy', ['Sz_total', 'which', 'L'])
plt.figure()
pyalps.plot.plot(sector_energies)
plt.xlabel('$J_1/J$')
plt.ylabel('$E_0$')
plt.legend(prop={'size':8})

# for each value of J1, L, we need to calculate the singlet and triplet gap:
# singlet: gap from abs. g.s. to first excited in S=0 sector
# triplet: gap from abs. g.s. to lowest in S=1 sector
grouped = pyalps.groupSets( pyalps.groupSets(pyalps.flatten(data), ['J1', 'L']), ['Sz_total'])

gaps = []
for J1g in grouped:
    totalmin = 1000
    for q in flatten(J1g):
        totalmin = min(totalmin, np.min(q.y))
    
    for Szg in J1g:
        allE = []
        for q in Szg:
            allE += list(q.y)
        allE = np.sort(allE)
        d = pyalps.DataSet()
        d.props = pyalps.dict_intersect([q.props for q in Szg])
        d.props['observable'] = 'gap'
        print(totalmin,d.props['Sz_total'])
        if d.props['Sz_total'] == 0:
            d.y = np.array([allE[1]-totalmin])
        else:
            d.y = np.array([allE[0]-totalmin])
        d.x = np.array([0])
        d.props['line'] = '.-'
        gaps.append(d)

gaps = pyalps.collectXY(gaps, 'J1', 'gap', ['Sz_total', 'L'])

plt.figure()
pyalps.plot.plot(gaps)
plt.xlabel('$J_1/J$')
plt.ylabel('$\Delta$')
plt.legend(prop={'size':8})

plt.show()
