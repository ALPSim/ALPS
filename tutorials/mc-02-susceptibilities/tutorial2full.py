# Please run all four other tutorials before running this one. 
# This tutorial relies on the results created in those tutorials

import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot

# load all files
data = pyalps.loadMeasurements(pyalps.getResultFiles(),'Susceptibility')

#flatten the hierarchical structure
data = pyalps.flatten(data)

# collect the susceptibility
susceptibility = pyalps.collectXY(data,x='T',y='Susceptibility',foreach=['MODEL','LATTICE'])

# assign labels to the data depending on the properties
for s in susceptibility:
  print s.props
  if s.props['LATTICE']=='chain lattice':
    s.props['label'] = "chain"
  elif s.props['LATTICE']=='ladder':
    s.props['label'] = "ladder"
  if s.props['MODEL']=='spin':
    s.props['label'] = "quantum " + s.props['label']
  elif s.props['MODEL']=='Heisenberg':
    s.props['label'] = "classical " + s.props['label']

#make plot
plt.figure()
pyalps.pyplot.plot(susceptibility)
plt.xlabel('Temperature $T/J$')
plt.ylabel('Susceptibility $\chi J$')
plt.ylim(0,0.25)
plt.legend()
plt.show()

