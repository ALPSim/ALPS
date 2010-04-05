import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot
import numpy as np

#prepare the input parameters
parms = [{ 
          'LATTICE'                   : "V15", 
          'LATTICE_LIBRARY'           : "v15-graph.xml", 
          'MODEL'                     : "spin",
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'local_S'                   : 0.5,
          'J'                         : 1
        }]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm8d',parms)
res = pyalps.runApplication('fulldiag',input_file)

#run the evaluation and load all the plots
data = pyalps.evaluateFulldiagVersusT(pyalps.getResultFiles(prefix='parm8d'),DELTA_T=0.05, T_MIN=0.05, T_MAX=5.0)

#make plot
for s in pyalps.flatten(data):
  plt.figure()
  plt.title("V$_{15}$ molecular magnet")
  pyalps.pyplot.plot(s)

