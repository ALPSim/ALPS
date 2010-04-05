import pyalps
import numpy as np
import matplotlib.pyplot as plt
import pyalps.pyplot

#prepare the input parameters
parms = { 
          'LATTICE'          : "chain lattice", 
          'LATTICE_LIBRARY'   : "lattices_dmrg.xml", 
          'MODEL'             : "spin",
          'L'                 : 10,
          't'                 : 1,
          'V'                 : 0,
          'SWEEPS'            : 10,
          'WAVEFUNCTION_FILE' : "psi.dat",
          'OUTPUT_LEVEL'      : 1
        }

#write the input file and run the simulation

input_file = pyalps.writeParameterFile('parm9a',parms)
res = pyalps.runApplication('simple_dmrg',input_file)

# load the text file and plot the wave function
raw = np.loadtxt(parms['WAVEFUNCTION_FILE']).transpose()
data = pyalps.DataSet()
data.x = raw[0]
data.y = raw[1]
data.props['xlabel'] = 'x'

plt.figure()
pyalps.pyplot.plot(data)
