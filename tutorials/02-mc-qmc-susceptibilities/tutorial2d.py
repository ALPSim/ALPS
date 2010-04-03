import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot

#prepare the input parameters
parms = []
for t in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0]:
    parms.append(
        { 
          'LATTICE'        : "ladder", 
          'MODEL'          : "spin",
          'local_S'        : 0.5,
          'T'              : t,
          'J0'             : 1 ,
          'J1'             : 1 ,
          'THERMALIZATION' : 10000,
          'SWEEPS'         : 150000,
          'L'              : 60,
          'ALGORITHM'      : "loop"
        }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm2d',parms)
res = pyalps.runApplication('loop',input_file)
output_file = res[1]
pyalps.evaluateLoop(output_file)

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm2d'),'Susceptibility')
susceptibility = pyalps.collectXY(data,x='T',y='Susceptibility')

#make plot
plt.figure()
pyalps.pyplot.plot(susceptibility)
plt.xlabel('Temperature $T/J$')
plt.ylabel('Susceptibility $\chi J$')
plt.ylim(0,0.22)
plt.title('Quantum Heisenberg ladder')
