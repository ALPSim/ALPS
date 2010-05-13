import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot

#prepare the input parameters
parms = []
for l in [4,8,16]: #[4,8,16]:
    for t in [2.8,2.7]:#[5.0,4.5,4.0,3.5,3.0,2.9,2.8,2.7]:
        parms.append(
            { 
              'LATTICE'        : "square lattice", 
              'T'              : t,
              'J'              : 1 ,
              'THERMALIZATION' : 1000,
              'SWEEPS'         : 40000, #400000,
              'UPDATE'         : "cluster",
              'MODEL'          : "Ising",
              'L'              : l
            }
    )
    for t in [2.6, 2.5, 2.4]: #[2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.2]:
        parms.append(
            { 
              'LATTICE'        : "square lattice", 
              'T'              : t,
              'J'              : 1 ,
              'THERMALIZATION' : 1000,
              'SWEEPS'         : 4000, #40000,
              'UPDATE'         : "cluster",
              'MODEL'          : "Ising",
              'L'              : l
            }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm7a',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7a'))

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7a'),['|Magnetization|', 'Connected Susceptibility', 'Specific Heat', 'Binder Cumulant', 'Binder Cumulant U2'])
magnetization_abs = pyalps.collectXY(data,x='T',y='|Magnetization|',foreach=['L'])
connected_susc = pyalps.collectXY(data,x='T',y='Connected Susceptibility',foreach=['L'])
spec_heat = pyalps.collectXY(data,x='T',y='Specific Heat',foreach=['L'])
binder_u4 = pyalps.collectXY(data,x='T',y='Binder Cumulant',foreach=['L'])
binder_u2 = pyalps.collectXY(data,x='T',y='Binder Cumulant U2',foreach=['L'])

#make plot
plt.figure()
pyalps.pyplot.plot(magnetization_abs)
plt.xlabel('Temperature $T$')
plt.ylabel('Magnetization $|m|$')
plt.title('2D Ising model')

plt.figure()
pyalps.pyplot.plot(connected_susc)
plt.xlabel('Temperature $T$')
plt.ylabel('Connected Susceptibility $\chi_c$')
plt.title('2D Ising model')

plt.figure()
pyalps.pyplot.plot(spec_heat)
plt.xlabel('Temperature $T$')
plt.ylabel('Specific Heat $c_v$')
plt.title('2D Ising model')

plt.figure()
pyalps.pyplot.plot(binder_u4)
plt.xlabel('Temperature $T$')
plt.ylabel('Binder Cumulant U4 $g$')
plt.title('2D Ising model')

plt.figure()
pyalps.pyplot.plot(binder_u2)
plt.xlabel('Temperature $T$')
plt.ylabel('Binder Cumulant U2 $g$')
plt.title('2D Ising model')
plt.show()

