import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot
import numpy as np

#prepare the input parameters
parms = [{ 
          'LATTICE'                   : "open ladder", 
          'MODEL'                     : "spin",
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'SWEEPS'                    : 10,
          'MAXSTATES'                 : 100,
          'L'                         : 10,
          'local_S'                   : 0.5,
          'Sz_total'                  : 0,
          'J'                         : 1,
          'h'                         : 0,
          'MEASURE_AVERAGE[Magnetization]'                      : 'Sz',
          'MEASURE_AVERAGE[Exchange]'                           : 'exchange',
          'MEASURE_LOCAL[Local magnetization]'                  : 'Sz',
          'MEASURE_CORRELATIONS[Diagonal spin correlations]='   : 'Sz',
          'MEASURE_CORRELATIONS[Offdiagonal spin correlations]' : 'Splus:Sminus'

        }]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm10d',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm10d'))

# print properties of the ground state:
for s in data[0]:
  if pyalps.size(s.y[0])==1:
    print s.props['observable'], ' : ', s.y[0]
  else:
    for (x,y) in zip(s.x,s.y[0]):
      print  s.props['observable'], '(', x, ') : ', y
