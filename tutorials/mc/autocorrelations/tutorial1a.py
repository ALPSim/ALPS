import pyalps
import matplotlib.pyplot as plt
import pyalps.pyplot

parms = []
for l in [2,4,8,16,32,48]:
    parms.append(
        { 
          'LATTICE'        : "square lattice", 
          'T'              : 2.269186,
          'J'              : 1 ,
          'THERMALIZATION' : 10000,
          'SWEEPS'         : 50000,
          'UPDATE'         : "local",
          'MODEL'          : "Ising",
          'L'              : l
        }
    )

input_file = pyalps.writeInputFiles('parm1a',parms)
#pyalps.execute('spinmc',input_file,Tmin=5)

binning = pyalps.loadBinningAnalysis(pyalps.getResultFiles(),'|Magnetization|')
binning = pyalps.flatten(binning)

for dataset in binning:
    dataset.props['label'] = 'L='+str(dataset.props['L'])

plt.figure()
plt.xlabel('binning level')
plt.ylabel('Error of |Magnetization|')
pyalps.pyplot.plot(binning)
plt.legend()

for dataset in binning:
    plt.figure()
    plt.title('Binning analysis for L='+str(dataset.props['L']))
    plt.xlabel('binning level')
    plt.ylabel('Error of |Magnetization|')
    pyalps.pyplot.plot(dataset)

# pyalps.writeGraceFile(binning,'parm1a.xmgr')



