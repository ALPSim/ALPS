# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 1994-2010 by Bela Bauer <bauerb@phys.ethz.ch>
#                            Ping Nang Ma
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the “Software”),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
# 
# ****************************************************************************

from __future__ import absolute_import
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
