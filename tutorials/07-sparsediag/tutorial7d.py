import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.pyplot

#prepare the input parameters
parms = []
for l in [4, 6, 8, 10]:
  for s in [0.5, 1]:
    for sz in [0, 1]:
      parms.append(
        { 
          'LATTICE'                   : "chain lattice", 
          'MODEL'                     : "spin",
          'local_S'                   : s,
          'J'                         : 1,
          'L'                         : l,
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'Sz_total'                  : sz
        }
      )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm7d',parms)
res = pyalps.runApplication('sparsediag',input_file)

#load all measurements for all states
data = pyalps.loadSpectra(pyalps.getResultFiles(prefix='parm7d'))

lengths = []
min_energies = {}

# extract the ground state energies over all momenta for every simulation
for sim in data:
  l = int(sim[0].props['L'])
  if l not in lengths: lengths.append(l)
  sz = int(sim[0].props['Sz_total'])
  s = float(sim[0].props['local_S'])
  all_energies = []
  for sec in sim:
    all_energies += list(sec.y)
  min_energies[(l,s,sz)]= np.min(all_energies)
  

# make a plot of the triplet gap as function of system size 
plt.figure()
for s in [0.5,1]:
  gapplot = pyalps.DataSet()
  gapplot.x = 1./np.sort(lengths)
  gapplot.y = [min_energies[(l,s,1)] -min_energies[(l,s,0)] for l in np.sort(lengths)]  
  gapplot.props['xlabel']='$1/L$'
  gapplot.props['ylabel']='Triplet gap $\Delta/J$'
  gapplot.props['label']='S='+str(s)
  pyalps.pyplot.plot(gapplot)

plt.legend()
plt.xlim(0,0.25)
plt.ylim(0,1.0)
plt.show()
