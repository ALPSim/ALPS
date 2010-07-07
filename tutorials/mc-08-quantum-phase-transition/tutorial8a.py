import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot
import numpy as np
import pyalps.fit_wrapper as fw
from math import sqrt

#prepare the input parameters
parms = []
for j2 in [0.,1.]:
    for t in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0]:
        parms.append(
            { 
              'LATTICE'        : "coupled ladders", 
              'LATTICE_LIBRARY': 'lattices.xml',
              'MODEL_LIBRARY'  : 'models.xml',
              'local_S'        : 0.5,
              'ALGORITHM'      : 'loop',
              'SEED'           : 0,
              'T'              : t,
              'J0'             : 1 ,
              'J1'             : 1,
              'J2'             : j2,
              'THERMALIZATION' : 5000,
              'SWEEPS'         : 50000, 
              'MODEL'          : "spin",
              'L'              : 8,
              'W'              : 4
            }
    )
    
#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm8a',parms)
res = pyalps.runApplication('loop',input_file)
output_file = res[1]
pyalps.evaluateLoop(output_file)

data = pyalps.loadMeasurements(pyalps.getResultFiles(pattern='parm8a.task*.out.h5'),['Staggered Susceptibility','Susceptibility'])
susc1=pyalps.collectXY(data,x='T',y='Susceptibility', foreach=['J2'])

lines = []
for data in susc1:
    pars = [fw.Parameter(1), fw.Parameter(1)]
    data.y= data.y[data.x < 1]
    data.x= data.x[data.x < 1]
    f = lambda self, x, pars: (pars[0]()/np.sqrt(x))*np.exp(-pars[1]()/x)
    fw.fit(None, f, pars, [v.mean for v in data.y], data.x)
    prefactor = pars[0].get()
    gap = pars[1].get()
    print prefactor,gap
    
    lines += plt.plot(data.x, f(None, data.x, pars))
    lines[-1].set_label('$J_2=%.4s$: $\chi = \\frac{%.4s}{T}\exp(\\frac{-%.4s}{T})$' % (data.props['J2'], prefactor,gap))

plt.figure()
pyalps.pyplot.plot(susc1)
plt.xlabel(r'$T$')
plt.ylabel(r'$\chi$')
plt.title('gap is %.4s' % gap)
plt.legend()
plt.show()
