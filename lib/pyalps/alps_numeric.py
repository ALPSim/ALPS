import numpy
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
  if (isinstance(obj,numpy.ndarray)) : \n\
    return numpy.OPERATION(obj)\n\
"


for operation in ["sq","sqrt","cb","cbrt","exp","log","sin","cos","tan","asin","acos","atan","sinh","cosh","tanh","asinh","acosh","atanh"]:
  function = global_function.replace("OPERATION",operation)
  exec function


  
  
