import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot

#prepare the input parameters
parms = []
for l in [32,48,64]:
    for t in [2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.30, 2.31, 2.32, 2.33, 2.34, 2.35]:
        parms.append(
            { 
              'LATTICE'        : "square lattice", 
              'T'              : t,
              'J'              : 1 ,
              'THERMALIZATION' : 5000,
              'SWEEPS'         : 150000,
              'UPDATE'         : "cluster",
              'MODEL'          : "Ising",
              'L'              : l
            }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm7b',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7b'))

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7b'),['|Magnetization|', 'Connected Susceptibility', 'Specific Heat', 'Binder Cumulant', 'Binder Cumulant U2'])
magnetization_abs = pyalps.collectXY(data,x='T',y='|Magnetization|',foreach=['L'])
connected_susc = pyalps.collectXY(data,x='T',y='Connected Susceptibility',foreach=['L'])
spec_heat = pyalps.collectXY(data,x='T',y='Specific Heat',foreach=['L'])
binder_u4 = pyalps.collectXY(data,x='T',y='Binder Cumulant',foreach=['L'])
binder_u2 = pyalps.collectXY(data,x='T',y='Binder Cumulant U2',foreach=['L'])

Tc=2.269
a=1

for d in binder_u4:
    d.x -= Tc
    d.x = d.x/Tc
    l = d.props['L']
    d.x = d.x * pow(float(l),a)
    
#find maximum in spec_heat and connected susc:
peaks=[]      
for d in connected_susc:
    peaks.append(np.max(d.y)) #can i somehow store the corresponding properties as well?


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

