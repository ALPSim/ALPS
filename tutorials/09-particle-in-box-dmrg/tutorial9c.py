import pyalps

#prepare the input parameters
parms = { 
          'LATTICE'          : "next-nearest chain lattice", 
          'LATTICE_LIBRARY'   : "lattices_dmrg.xml", 
          'MODEL'             : "spin",
          'L'                 : 20,
          't'                 : 1,
          't1'                : 0.01,
          'K'                 : '2*3.1415927*2/L',
          'V'                 : 'cos(K*x)+cos(3*K*x)/9+cos(5*K*x)/25+cos(7*K*x)/36-2',
          'SWEEPS'            : 100,
          'WAVEFUNCTION_FILE' : "psi2.dat",
          'OUTPUT_LEVEL'      : 1
        }

#write the input file and run the simulation

input_file = pyalps.writeParameterFile('parm9c',parms)
res = pyalps.runApplication('simple_dmrg',input_file)

# load the text file and plot the wave function
raw = np.loadtxt(parms['WAVEFUNCTION_FILE']).transpose()
data = pyalps.DataSet()
data.x = raw[0]
data.y = raw[1]
data.props['xlabel'] = 'x'

plt.figure()
pyalps.pyplot.plot(data)
