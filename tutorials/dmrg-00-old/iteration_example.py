import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.pyplot

#load all measurements for all states
data = pyalps.loadIterationMeasurements(pyalps.getResultFiles(prefix='parm10a'))
energy_iterations = pyalps.collectXY(data,x='iteration',y='Energy')

#make plot
plt.figure()
energy_iterations[0].props['line']='+-'
pyalps.pyplot.plot(energy_iterations)
plt.xlabel('Iteration')
plt.ylabel('Energy $E$')
plt.title('Energy Convergence in DMRG')
plt.show()



