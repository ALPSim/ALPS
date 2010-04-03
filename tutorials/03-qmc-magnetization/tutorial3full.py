# Please run the two other tutorials before running this one. 
# This tutorial relies on the results created in those tutorials

import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot

# load all files
data = pyalps.loadMeasurements(pyalps.getResultFiles(),'Magnetization Density')

#flatten the hierarchical structure
data = pyalps.flatten(data)

#load the magnetization and collect it as function of field h
magnetization = pyalps.collectXY(data,x='h',y='Magnetization Density',foreach=['LATTICE'])

#make plot
plt.figure()
pyalps.pyplot.plot(magnetization)
plt.xlabel('Field $h$')
plt.ylabel('Magnetization $m$')
plt.ylim(0.0,0.5)
plt.legend()
