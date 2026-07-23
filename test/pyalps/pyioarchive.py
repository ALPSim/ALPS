 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 #                                                                                 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   #
 #                                                                                 #
 # SPDX-License-Identifier: MIT
 #                                                                                 #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import pyalps.hdf5 as h5

o = h5.archive('blubb', 'w')
o['/a'] = 0
del o

i = h5.archive('blubb', 'r')
o = h5.archive('blubb', 'w')
o['/a'] = 0
del o

del i
o = h5.archive('blubb', 'w')
o['/a'] = 0
del o