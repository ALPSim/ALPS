import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot

#prepare the input parameters
parms = [{ 
          'LATTICE'        : "simple cubic lattice", 
          'MODEL'          : "spin",
          'local_S'        : 0.5,
          'L'              : 4,
          'J'              : 1 ,
          'CUTOFF'         : 500
        }]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm6c',parms)
res = pyalps.runApplication('qwl',input_file)

#run the evaluation and load all the plots
data = pyalps.evaluateQWL(pyalps.getResultFiles(prefix='parm6c'),DELTA_T=0.05, T_MIN=0.5, T_MAX=5.0)

#make plot
for s in pyalps.flatten(data):
  plt.figure()
  plt.title("Cubic lattice Heisenberg antiferromagnet L=4")
  pyalps.pyplot.plot(s)
