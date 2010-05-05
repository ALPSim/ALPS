import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot

#prepare the input parameters
parms = []
for h in [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.5]:
    parms.append(
        { 
          'LATTICE'        : "chain lattice", 
          'MODEL'          : "spin",
          'local_S'        : 0.5,
          'T'              : 0.08,
          'J'              : 1 ,
          'THERMALIZATION' : 1000,
          'SWEEPS'         : 10000,
          'L'              : 20,
          'h'              : h
        }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm3a',parms)
res = pyalps.runApplication('dirloop_sse',input_file,Tmin=5)

#load the magnetization and collect it as function of field h
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm3a'),'Magnetization Density')
magnetization = pyalps.collectXY(data,x='h',y='Magnetization Density')

#make plot
plt.figure()
pyalps.pyplot.plot(magnetization)
plt.xlabel('Field $h$')
plt.ylabel('Magnetization $m$')
plt.ylim(0.0,0.5)
plt.title('Quantum Heisenberg chain')
plt.show()
