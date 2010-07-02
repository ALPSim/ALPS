import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot
import numpy as np

#prepare the input parameters
parms = []
for l in [8,10,12,16]:
    for j2 in [1.8,1.85,1.9,1.95,2.,2.05,2.1]:
        parms.append(
            { 
              'LATTICE'        : "coupled ladders", 
              'LATTICE_LIBRARY': 'lattices.xml',
              'MODEL_LIBRARY'  : 'models.xml',
              'local_S'        : 0.5,
              'ALGORITHM'      : 'loop',
              'SEED'           : 0,
              'BETA'           : 2*l,
              'J0'             : 1 ,
              'J1'             : 1,
              'J2'             : j2,
              'THERMALIZATION' : 5000,
              'SWEEPS'         : 50000, 
              'MODEL'          : "spin",
              'L'              : l,
              'W'              : l/2
            }
    )
    
#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm8b',parms)
res = pyalps.runApplication('loop',input_file)
output_file = res[1]
pyalps.evaluateLoop(output_file)

data = pyalps.loadMeasurements(pyalps.getResultFiles(pattern='parm8a.task*.out.h5'),['Binder Ratio of Staggered Magnetization','Stiffness'])

binder=pyalps.collectXY(data,x='J2',y='Binder Ratio of Staggered Magnetization', foreach=['L'])
stiffness =pyalps.collectXY(data,x='J2',y='Stiffness', foreach=['L'])

for q in stiffness:
    q.y = q.y*q.props['L']

#make plot    
plt.figure()
pyalps.pyplot.plot(stiffness)
plt.xlabel(r'$J2$')
plt.ylabel(r'Stiffness $\rho_s L$')
plt.title('coupled ladders')

plt.figure()
pyalps.pyplot.plot(binder)
plt.xlabel(r'$J_2$')
plt.ylabel(r'$g(m_s)$')
plt.title('coupled ladders')
plt.show()
