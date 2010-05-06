import pyalps
import pyalps.pyplot
import numpy as np
import matplotlib.pyplot as plt
import copy
import math

# Some general parameters
parms_ = {
	'LATTICE'              : "chain lattice",
	'MODEL'                : "spin",
	'local_S'              : 0.5,
	'J'                    : 1,
	'NUMBER_EIGENVALUES'   : 5,
	'CONSERVED_QUANTUMNUMBER' : 'Sz',
	'Sz_total' : 0
}

prefix = 'heisenberg'
parms = []
for L in [10,12]:
	parms_.update({'L':L})
	parms.append(copy.deepcopy(parms_))

input_file = pyalps.writeInputFiles(prefix,parms)
res = pyalps.runApplication('sparsediag', input_file)

data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix=prefix))

E0 = {}
E1 = {}
for Lsets in data:
	L = pyalps.flatten(Lsets)[0].props['L']
	# Make a big list of all energy values
	allE = []
	for q in pyalps.flatten(Lsets):
		allE += list(q.y)
	allE = np.sort(allE)
	E0[L] = allE[0]
	E1[L] = allE[1]

for q in pyalps.flatten(data):
	L = q.props['L']
	q.y = (q.y-E0[L])/(E1[L]-E0[L]) * (1./2.)

spectrum = pyalps.collectXY(data, 'TOTAL_MOMENTUM', 'Energy', foreach=['L'])

for SD in [0.5, 1, 1+0.5, 1+1]:
	d = pyalps.DataSet()
	d.x = np.array([0,4])
	d.y = SD+0*d.x
	# d.props['label'] = str(SD)
	spectrum += [d]

pyalps.pyplot.plot(spectrum)

plt.legend(prop={'size':8})
plt.xlabel("$k$")
plt.ylabel("E_0")

plt.xlim(-0.02, math.pi+0.02)

plt.show()
