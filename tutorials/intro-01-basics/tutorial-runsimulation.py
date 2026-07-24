from __future__ import print_function
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

# set the input file 
input_file = 'parm1.in.xml'

# and run the simulation
result_files = pyalps.runApplication('spinmc',input_file,Tmin=5,writexml=True)
print("The return value was", result_files[0])
print("The results are stored in the file", result_files[1])
