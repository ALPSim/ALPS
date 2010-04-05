import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot
import numpy as np

#prepare the input parameters
parms = [{ 
          'LATTICE'                   : "open chain lattice", 
          'MODEL'                     : "hardcore boson",
          'CONSERVED_QUANTUMNUMBERS'  : 'N',
          'SWEEPS'                    : 5,
          'STATES'                    : "20,20,40,40,60,60,80,80,100,100",
          'L'                         : 10,
          'N_total'                   : 2,
          't'                         : 1,
          'V'                         : 1.3
        }]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm10e',parms)
res = pyalps.runApplication('dmrg',input_file,writexml=True)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm10e'))

# print properties of the ground state:
for s in data[0]:
  if pyalps.size(s.y[0])==1:
    print s.props['observable'], ' : ', s.y[0]
  else:
    for (x,y) in zip(s.x,s.y[0]):
      print  s.props['observable'], '(', x, ') : ', y
