import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot
import numpy as np

#prepare the input parameters
parms = []
for h in [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]:
    parms.append(
        { 
          'LATTICE'                        : "open chain lattice", 
          'MODEL'                          : "spin",
          'CONSERVED_QUANTUMNUMBERS'       : 'Sz',
          'SWEEPS'                         : 4,
          'MAXSTATES'                      : 100,
          'L'                              : 10,
          'local_S'                        : 0.5,
          'J'                              : 1,
          'h'                              : h,
          'MEASURE_AVERAGE[Magnetization]' : 'Sz'
        }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm10b',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load the magnetization and collect it as function of field h
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm10b'))
magnetization = pyalps.collectXY(data,x='h',y='Magnetization')

magnetization[0].y /= magnetization[0].props['L']
#make plot
plt.figure()
pyalps.pyplot.plot(magnetization)
plt.xlabel('Field $h$')
plt.ylabel('Magnetization density $m$')
plt.ylim(0.0,0.5)
plt.title('Quantum Heisenberg chain')
