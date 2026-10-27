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

# Python math and NumPy spell some of these operations differently, have no
# sq/cb at all, and math.cbrt needs Python 3.11. Other types (alea
# observables) provide every operation as a method.
_scalar_ops = {"sq": lambda x: x*x, "cb": lambda x: x*x*x, "cbrt": lambda x: pm.copysign(abs(x)**(1./3.), x)}
_array_ops = {"sq": np.square, "cb": lambda x: x*x*x,
              "asin": np.arcsin, "acos": np.arccos, "atan": np.arctan,
              "asinh": np.arcsinh, "acosh": np.arccosh, "atanh": np.arctanh}

def _make_dispatcher(operation):
  scalar_op = _scalar_ops.get(operation) or getattr(pm, operation)
  array_op = _array_ops.get(operation) or getattr(np, operation)
  def dispatcher(obj):
    # math library
    if (isinstance(obj,(float,int))):
      return scalar_op(obj)
    # numpy array
    if (isinstance(obj,np.ndarray)):
      return array_op(obj)
    # other types
    return getattr(obj, operation)()
  dispatcher.__name__ = operation
  return dispatcher

for operation in ["sq", "sqrt", "cb", "cbrt", "exp", "log", "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh", "asinh", "acosh", "atanh"]:
  globals()[operation] = _make_dispatcher(operation)
