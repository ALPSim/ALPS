# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 1994-2010 by Bela Bauer <bauerb@phys.ethz.ch>
#                            Ping Nang Ma
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import numpy as np
from pyalps.alea import *
import math as pm

global_function = "def OPERATION(obj): \n\
  # math library \n\
  if (isinstance(obj,float)): \n\
    return pm.OPERATION(obj) \n\
  if (isinstance(obj,int)): \n\
    return pm.OPERATION(obj) \n\
\n\
  # numpy array \n\
  if (isinstance(obj,np.ndarray)) : \n\
    return np.OPERATION(obj)\n\
  # other types \n\
  return obj.OPERATION() \n\
\n\
"

for operation in ["sq", "sqrt", "cb", "cbrt", "exp", "log", "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh", "asinh", "acosh", "atanh"]:
  function = global_function.replace("OPERATION",operation)
  exec(function)
