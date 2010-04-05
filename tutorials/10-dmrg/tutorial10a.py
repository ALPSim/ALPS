import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot
import numpy as np

#prepare the input parameters
parms = [{ 
          'LATTICE'                   : "open chain lattice", 
          'MODEL'                     : "spin",
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'SWEEPS'                    : 4,
          'MAXSTATES'                 : 100,
          'L'                         : 10,
          'local_S'                   : 0.5,
          'Sz_total'                  : 0,
          'J'                         : 1,
          'h'                         : 3,
          'NUMBER_EIGENVALUES'        : 2,
          'MEASURE_AVERAGE[Magnetization]'                      : 'Sz',
          'MEASURE_AVERAGE[Exchange]'                           : 'exchange',
          'MEASURE_LOCAL[Local magnetization]'                  : 'Sz',
          'MEASURE_CORRELATIONS[Diagonal spin correlations]='   : 'Sz',
          'MEASURE_CORRELATIONS[Offdiagonal spin correlations]' : 'Splus:Sminus'

        }]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm10a',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm10a'))

# print properties of all eigenvectors:
for index in range(0,2):
  for s in data[0]:
    if index < len(s.y):
      if pyalps.size(s.y[index])==1:
        print s.props['observable'], ' : ', s.y[index]
      else:
        for (x,y) in zip(s.x,s.y[index]):
          print  s.props['observable'], '(', x, ') : ', y
