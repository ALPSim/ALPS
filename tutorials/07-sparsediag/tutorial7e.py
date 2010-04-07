import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.pyplot

#prepare the input parameters

parms=[]
for sz in [0, 1]:
      parms.append(
        {
          'LATTICE'                   : "ladder", 
          'MODEL'                     : "spin",
          'local_S'                   : 0.5,
          'J0'                        : 1,
          'J1'                        : 1,
          'L'                         : 8,
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'Sz_total'                  : sz,
          'NUMBER_EIGENVALUES'        : 10
        }
      )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm7e',parms)
res = pyalps.runApplication('sparsediag',input_file)

#load all measurements for all states
data = pyalps.loadSpectra(pyalps.getResultFiles(prefix='parm7e'))

energies=[]
# get ground state energy
for s in pyalps.flatten(data):
    energies += list(s.y)

groundstate_energy = np.min(energies)

#plot spectra
spectrumplot = []
for set in data:
  plot = pyalps.DataSet()
  plot.props['label']='Sz='+str(set[0].props['Sz_total'])
  plot.props['line'] = 'scatter'
  for s in set:
    plot.x = np.concatenate((plot.x,np.array([s.props['TOTAL_MOMENTUM'] for i in range(0,len(s.y))])))
    plot.y = np.concatenate((plot.y,s.y - groundstate_energy))
  spectrumplot.append(plot)
  
#spectrumplot[0] = pyalps.subtract_spectrum(spectrumplot[0],spectrumplot[1],tolerance=1e-12)

plt.figure()
pyalps.pyplot.plot(spectrumplot)
plt.legend()
plt.xlim(0,2*3.1416)
plt.ylabel('Energy')
plt.ylim(0,3)




