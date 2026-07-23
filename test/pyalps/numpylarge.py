from __future__ import print_function
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

import numpy as np
import pyalps.hdf5 as hdf5
import os

ar=hdf5.archive('foo%d.h5', 'al')
s=2**10

while s < 2**29:
    print(s)
    a = np.empty(s)
    ar[str(s)] = a
    s *= 2

i = 0
while os.path.isfile('foo%d.h5'%i):
    os.remove('foo%d.h5'%i)
    i += 1
