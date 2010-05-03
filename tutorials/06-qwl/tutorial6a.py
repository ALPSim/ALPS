import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot

#prepare the input parameters
parms = [{ 
          'LATTICE'        : "chain lattice", 
          'MODEL'          : "spin",
          'local_S'        : 0.5,
          'L'              : 40,
          'J'              : -1 ,
          'CUTOFF'         : 1000
        }]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm6a',parms)
res = pyalps.runApplication('qwl',input_file)

#run the evaluation and load all the plots
data = pyalps.evaluateQWL(pyalps.getResultFiles(prefix='parm6a'),DELTA_T=0.1, T_MIN=0.1, T_MAX=10.0)

#make plot
for s in pyalps.flatten(data):
  plt.figure()
  plt.title("Ferromagnetic Heisenberg chain")
  pyalps.pyplot.plot(s)
plt.show()
