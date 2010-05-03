import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot
import numpy as np

#prepare the input parameters
parms = [{ 
          'LATTICE'                   : "chain lattice", 
          'MODEL'                     : "spin",
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'local_S'                   : 1,
          'J'                         : 1,
          'L'                         : 8
        }]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm8a',parms)
res = pyalps.runApplication('fulldiag',input_file)

#run the evaluation and load all the plots
data = pyalps.evaluateFulldiagVersusT(pyalps.getResultFiles(prefix='parm8a'),DELTA_T=0.1, T_MIN=0.1, T_MAX=10.0)

#make plot
for s in pyalps.flatten(data):
  plt.figure()
  plt.title("Antiferromagnetic Heisenberg chain")
  pyalps.pyplot.plot(s)


# make a plot of the spectrum, first load all measurements for all states
spectrum = pyalps.loadSpectra(pyalps.getResultFiles(prefix='parm8a'))

energies=[]
# get ground state energy
for set in spectrum:
  for s in set:
    energies += list(s.y)

groundstate_energy = np.min(energies)

#collect spctra for each Sz into one DataSet in a Python dictionary
spectrumplot = {}
for sz in [0,1,2]:
    spectrumplot[sz] = pyalps.DataSet()
    spectrumplot[sz].props['label']='Sz='+str(sz)
    spectrumplot[sz].props['line'] = 'scatter'

for s in spectrum[0]:
  sz = int(s.props['Sz'])
  if sz in spectrumplot:
    spectrumplot[sz].x = np.concatenate((spectrumplot[sz].x,np.array([s.props['TOTAL_MOMENTUM'] for i in range(0,len(s.y))])))
    spectrumplot[sz].y = np.concatenate((spectrumplot[sz].y,s.y - groundstate_energy))

plt.figure()
pyalps.pyplot.plot(spectrumplot.values())
plt.legend()
plt.xlim(0,2*3.1416)
plt.ylabel('Energy')
plt.ylim(0,3)
plt.show()






