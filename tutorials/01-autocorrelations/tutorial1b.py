import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot

#prepare the input parameters
parms = []
for l in [2,4,8,16,32,48]:
    parms.append(
        { 
          'LATTICE'        : "square lattice", 
          'T'              : 2.269186,
          'J'              : 1 ,
          'THERMALIZATION' : 1000,
          'SWEEPS'         : 100000,
          'UPDATE'         : "cluster",
          'MODEL'          : "Ising",
          'L'              : l
        }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm1b',parms)
pyalps.execute('spinmc',input_file,Tmin=5)

#load the binning analysis for the absolute value of the magnetization
binning = pyalps.loadBinningAnalysis(pyalps.getResultFiles(prefix='parm1b'),'|Magnetization|')
binning = pyalps.flatten(binning)

#make one plot with all data
for dataset in binning:
    dataset.props['label'] = 'L='+str(dataset.props['L'])

plt.figure()
plt.title('Binning analysis for cluster updates')
plt.xlabel('binning level')
plt.ylabel('Error of |Magnetization|')
pyalps.pyplot.plot(binning)
plt.legend()


# make individual plots for each system size
for dataset in binning:
    plt.figure()
    plt.title('Binning analysis for L='+str(dataset.props['L']))
    plt.xlabel('binning level')
    plt.ylabel('Error of |Magnetization|')
    pyalps.pyplot.plot(dataset)



