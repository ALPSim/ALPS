import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot
import copy

#prepare the input parameters
parms=[]

for (l,c) in [(4,500),( 6,1000)]:
    parms.append(
       { 
          'LATTICE'        : "simple cubic lattice", 
          'MODEL'          : "spin",
          'local_S'        : 0.5,
          'L'              : l,
          'J'              : 1 ,
          'CUTOFF'         : c
        }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm6d',parms)
pyalps.runApplication('qwl',input_file)

#run the evaluation and load all the plots
results = pyalps.evaluateQWL(pyalps.getResultFiles(prefix='parm6d'),DELTA_T=0.05, T_MIN=0.5, T_MAX=1.5)

#extract just the staggered structure factor S(Q) and rescale it by L^{-2+\eta}
data = []
for s in pyalps.flatten(results):
  if s.props['ylabel']=='Staggered Structure Factor per Site':
    print 'yes'
    d = copy.deepcopy(s) # make a deep copy to not change the original
    l = s.props['L']
    d.props['label']='L='+str(l)
    d.y = d.y * pow(float(l),-1.97)
    data.append(d)

#make plot
plt.figure()
plt.title("Scaling plot for cubic lattice Heisenberg antiferromagnet")
pyalps.pyplot.plot(data)
plt.legend()
plt.xlabel('Temperature $T/J$')
plt.ylabel('$S(\pi,\pi,\pi) L^{-2+\eta}$')
