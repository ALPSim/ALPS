import pyalps
import pyalps.pyplot
from pyalps.dict_intersect import dict_intersect
import numpy as np
import matplotlib.pyplot as plt
import copy
import math

# Some general parameters
parms_ = {
	'LATTICE'              : "nnn chain lattice",
	'MODEL'                : "spin",
	'local_S'              : 0.5,
	'J'                    : 1,
	'NUMBER_EIGENVALUES'   : 2,
	'CONSERVED_QUANTUMNUMBER' : 'Sz'
}

prefix = 'alps-nnn-heisenberg'
parms = []
for L in [10]:
    for Szt in [0,1]:
    	for J1 in np.linspace(0,1,11):
    		parms_.update({'Sz_total':Szt, 'J1':J1, 'L':L})
    		parms.append(copy.deepcopy(parms_))

input_file = pyalps.writeInputFiles(prefix,parms)
# res = pyalps.runApplication('sparsediag', input_file)
res = pyalps.runApplication('sparsediag', input_file, MPI=4)
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix=prefix))

sector_E = []
grouped = pyalps.groupSets(pyalps.flatten(data), ['Sz_total', 'J1', 'L'])
for group in grouped:
    allE = []
    for q in group:
        allE += list(q.y)
    allE = np.sort(allE)
    
    d = pyalps.DataSet()
    d.props = dict_intersect([q.props for q in group])
    d.x = np.array([0])
    d.y = np.array([allE[0]])
    d.props['which'] = 'gs'
    sector_E.append(d)
    
    d2 = copy.deepcopy(d)
    d2.y = np.array([allE[1]])
    d2.props['which'] = 'fe'
    sector_E.append(d2)
    
    d3 = copy.deepcopy(d)
    d3.y = np.array([allE[2]])
    d3.props['which'] = 'se'
    sector_E.append(d3)

sector_energies = pyalps.collectXY(sector_E, 'J1', 'Energy', ['Sz_total', 'which', 'L'])
plt.figure()
pyalps.pyplot.plot(sector_energies)
plt.legend(prop={'size':8})

grouped = pyalps.groupSets( pyalps.groupSets(pyalps.flatten(data), ['Sz_total']), ['J1', 'L'])

gaps = []
for J1g in grouped:
    for Szg in J1g:
        allE = []
        for q in Szg:
            allE += list(q.y)
        allE = np.sort(allE)
        d = pyalps.DataSet()
        d.props = dict_intersect([q.props for q in Szg])
        d.props['observable'] = 'gap'
        # they may be degenerate
        if allE[1]-allE[0] > 1e-6:
            d.y = np.array([allE[1]-allE[0]])
        else:
            d.y = np.array([allE[2]-allE[0]])
        d.x = np.array([0])
        gaps.append(d)

nicer = pyalps.collectXY(gaps, 'J1', 'gap', ['Sz_total', 'L'])

plt.figure()
pyalps.pyplot.plot(nicer)
plt.legend(prop={'size':8})

