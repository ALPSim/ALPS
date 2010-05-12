import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.pyplot
import pyalps.fit_wrapper as fw

#prepare the input parameters
parms = []
for l in [4, 6, 8, 10, 12, 14]:
  for sz in [0, 1]:
      parms.append(
        { 
          'LATTICE'                   : "chain lattice", 
          'MODEL'                     : "spin",
          'local_S'                   : 1,
          'J'                         : 1,
          'L'                         : l,
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'Sz_total'                  : sz
        }
      )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm2a',parms)
res = pyalps.runApplication('sparsediag',input_file) # ,MPI=4)

#load all measurements for all states
data = pyalps.loadSpectra(pyalps.getResultFiles(prefix='parm2a'))

lengths = []
min_energies = {}

# extract the ground state energies over all momenta for every simulation
for sim in data:
  l = int(sim[0].props['L'])
  if l not in lengths: lengths.append(l)
  sz = int(sim[0].props['Sz_total'])
  all_energies = []
  for sec in sim:
    all_energies += list(sec.y)
  min_energies[(l,sz)]= np.min(all_energies)
  

# make a plot of the triplet gap as function of system size   
gapplot = pyalps.DataSet()
gapplot.x = 1./np.sort(lengths)
gapplot.y = [min_energies[(l,1)] -min_energies[(l,0)] for l in np.sort(lengths)]  
gapplot.props['xlabel']='$1/L$'
gapplot.props['ylabel']='Triplet gap $\Delta/J$'
gapplot.props['label']='S=1'
gapplot.props['line']='.'

plt.figure()
pyalps.pyplot.plot(gapplot)
plt.legend()
plt.xlim(0,0.25)
plt.ylim(0,1.0)

d0 = fw.Parameter(0.411)
L0 = fw.Parameter(1000)
a = fw.Parameter(1)
f = lambda self, x: d0()+a()*np.exp(-x/L0())
# we fit only a range from 8 to 14
fw.fit(None, f, (d0,L0,a), np.array(gapplot.y)[2:], np.sort(lengths)[2:])

x = np.linspace(0, 1./min(lengths), 100)
plt.plot(x, f(None, 1/x))

plt.show()


