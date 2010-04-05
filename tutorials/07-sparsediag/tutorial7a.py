import pyalps

#prepare the input parameters
parms = [{ 
          'LATTICE'                   : "chain lattice", 
          'MODEL'                     : "spin",
          'local_S'                   : 1,
          'J'                         : 1,
          'L'                         : 4,
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'MEASURE_LOCAL[Local magnetization]'                  : 'Sz',
          'MEASURE_STRUCTURE_FACTOR[Structure Factor S]'        : 'Sz',
          'MEASURE_CORRELATIONS[Diagonal spin correlations]='   : 'Sz',
          'MEASURE_CORRELATIONS[Offdiagonal spin correlations]' : 'Splus:Sminus'
        }]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm7a',parms)
res = pyalps.runApplication('sparsediag',input_file)

#load all measurements for all states
data = pyalps.loadEigenstateMeasurements(pyalps.getResultFiles(prefix='parm7a'))

# print properties of ground states in all sectors:
for sector in data[0]:
  print '\nSector with Sz =', sector[0].props['Sz'], 
  print 'and k =', sector[0].props['TOTAL_MOMENTUM']
  for s in sector:
    if pyalps.size(s.y[0])==1:
      print s.props['observable'], ' : ', s.y[0]
    else:
      for (x,y) in zip(s.x,s.y[0]):
        print  s.props['observable'], '(', x, ') : ', y

