import pyalps

#prepare the input parameters
parms = [{ 
          'LATTICE'                   : "square lattice", 
          'MODEL'                     : "spin",
          'MEASURE[Correlations]'     : True,
          'MEASURE[Structure Factor]' : True,
          'MEASURE[Green Function]'   : True,
          'local_S'                   : 0.5,
          'T'                         : 0.3,
          'J'                         : 1 ,
          'THERMALIZATION'            : 10000,
          'SWEEPS'                    : 500000,
          'L'                         : 4,
          'h'                         : 0.1
        }]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm4',parms)
res = pyalps.runApplication('dirloop_sse',input_file,Tmin=5)

#load the magnetization and collect it as function of field h
data = pyalps.loadMeasurements(pyalps.getResultFiles())

# print all measurements
for s in pyalps.flatten(data):
  if len(s.x)==1:
    print s.props['observable'], ' : ', s.y[0]
  else:
    for (x,y) in zip(s.x,s.y):
      print  s.props['observable'], x, ' : ', y
