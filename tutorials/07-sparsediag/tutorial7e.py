import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.pyplot

#prepare the input parameters
parms = [{ 
          'LATTICE'                   : "ladder", 
          'MODEL'                     : "spin",
          'local_S'                   : 0.5,
          'J0'                        : 1,
          'J1'                        : 1,
          'L'                         : 8,
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz'
        }]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm7e',parms)
res = pyalps.runApplication('sparsediag',input_file)

#load all measurements for all states
data = pyalps.loadSpectra(pyalps.getResultFiles(prefix='parm7e'))


plt.figure()
pyalps.pyplot.plot(gapplot)
plt.legend()
plt.xlim(0,0.25)
plt.ylim(0,1.0)




