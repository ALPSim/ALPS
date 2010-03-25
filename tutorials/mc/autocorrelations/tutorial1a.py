import pyalps

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
    
pyalps.writeInputFiles('parm1a',parms)

pyalps.execute('spinmc','parm1a.in.xml',Tmin=5)

binning = pyalps.loadBinningAnalysis(getResultFiles(),'|Magnetization|')



