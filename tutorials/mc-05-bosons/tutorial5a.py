import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot

#prepare the input parameters
parms = []
for t in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]:
    parms.append(
        { 
          'LATTICE'        : "square lattice", 
          'MODEL'          : "boson Hubbard",
          'T'              : 0.1,
          'L'              : 4 ,
          't'              : t ,
          'mu'             : 0.5,
          'U'              : 1.0 ,
          'NONLOCAL'       : 0 ,
          'Nmax'           : 2 ,
          'THERMALIZATION' : 10000,
          'SWEEPS'         : 500000
        }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm5a',parms)
res = pyalps.runApplication('worm',input_file,Tmin=5)

#load the magnetization and collect it as function of field h
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm5a'),'Stiffness')
rhos = pyalps.collectXY(data,x='t',y='Stiffness')

#make plot
plt.figure()
pyalps.pyplot.plot(rhos)
plt.xlabel('Hopping $t/U$')
plt.ylabel('Superfluid density $\\rho _s$')
plt.show()
