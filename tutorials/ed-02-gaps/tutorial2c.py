import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.pyplot

data = pyalps.loadSpectra(pyalps.getResultFiles(prefix='parm2a'))
data += pyalps.loadSpectra(pyalps.getResultFiles(prefix='parm2b'))

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
