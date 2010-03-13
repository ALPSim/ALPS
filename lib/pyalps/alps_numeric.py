# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 1994-2010 by Bela Bauer <bauerb@phys.ethz.ch>
#                            Ping Nang Ma
#
# This software is part of the ALPS libraries, published under the ALPS
# Library License; you can use, redistribute it and/or modify it under
# the terms of the license, either version 1 or (at your option) any later
# version.
#  
# You should have received a copy of the ALPS Library License along with
# the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
# available from http://alps.comp-phys.org/.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
# FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
# 
# ****************************************************************************


import numpy as np
import math
from pyalea import *


global_function = "def OPERATION(obj): \n\
  # rank-0 objects \n\
\n\
  # math library \n\
  if (isinstance(obj,float)): \n\
    return math.OPERATION(obj) \n\
  if (isinstance(obj,int)): \n\
    return math.OPERATION(obj) \n\
\n\
  # value_with_error library \n\
  if (isinstance(obj,value_with_error)): \n\
    return value_with_error.OPERATION(obj) \n\
  if (isinstance(obj,vector_with_error)): \n\
    return vector_with_error.OPERATION(obj) \n\
  if (isinstance(obj,vector_of_value_with_error)): \n\
    return vector_of_value_with_error.OPERATION(obj) \n\
\n\
  # numpy array \n\
  if (isinstance(obj,np.ndarray)) : \n\
    return np.OPERATION(obj)\n\
"


for operation in ["sq","sqrt","cb","cbrt","exp","log","sin","cos","tan","asin","acos","atan","sinh","cosh","tanh","asinh","acosh","atanh"]:
  function = global_function.replace("OPERATION",operation)
  exec function


  
  
