from __future__ import print_function
# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Lukas Gamper <gamperl@gmail.com>
#                       Matthias Troyer <troyer@itp.phys.ethz.ch>
#
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps as alps
import pyalps.hdf5 as h5
import pyalps.alea as alea

iar = h5.archive('loadobs.h5', 'r')

for name in iar.list_children('/simulation/results'):
    if iar.is_scalar('/simulation/results/' + alps.hdf5_name_encode(name) + '/mean/value'):
        obs = alea.MCScalarData()
    else:
        obs = alea.MCVectorData()
    obs.load('loadobs.h5', '/simulation/results/' + alps.hdf5_name_encode(name))
    print(name + ": " + str(obs))
