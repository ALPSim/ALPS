# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Matthias Troyer <troyer@phys.ethz.ch> 
# 
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps
import matplotlib.pyplot as plt
import pyalps.plot

#prepare the input parameters
parms = [{ 
          'LATTICE'                   : "double dimer", 
          'LATTICE_LIBRARY'           : "dd-graph.xml", 
          'MODEL'                     : "dimerized spin",
          'MODEL_LIBRARY'             : "model-dspin.xml",
          'CONSERVED_QUANTUMNUMBERS'  : 'Sz',
          'local_S0'                  : 1,
          'local_S1'                  : 0.5,
          'J0'                        : 1,
          'J1'                        : 0.4
        }]

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('ed06c',parms)
res = pyalps.runApplication('fulldiag',input_file)

#run the evaluation and load all the plots
data = pyalps.evaluateFulldiagVersusH(pyalps.getResultFiles(prefix='ed06c'),T = 0.02, DELTA_H=0.025, H_MIN=0., H_MAX=4.0)

#make plot
for s in pyalps.flatten(data):
  plt.figure()
  pyalps.plot.plot(s)

plt.show()
