#
# --> This is a float with error data type, written in python language.
#    -> It is basically an extension of the float data type, that takes care of uncorrelated error propagation.
#    -> There is currently no documentation on it at the moment.
#
# --> Written by Ping Nang (Tama) Ma, (pingnang@itp.phys.ethz.ch), ETH Zurich (Nov 2009)
#    
# 
#!/usr/bin/python

from sys  import stdin
from math import *

class FloatWithError:


  def __init__(self,mean_=0,error_=0):
    self.mean  = mean_
    self.error = error_

  def __str__(self):
    return str(self.mean) + '\t' + str(self.error)
  def __expr__(self):
    return expr(self.mean) + '\t' + expr(self.error)
  def __repr__(self):
	return self.__str__()
  
  def __add__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    return FloatWithError(x__.mean+y__.mean,sqrt(x__.error*x__.error+y__.error*y__.error))

  def __radd__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    return FloatWithError(y__.mean+x__.mean,sqrt(y__.error*y__.error+x__.error*x__.error))

  def __sub__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    return FloatWithError(x__.mean-y__.mean,sqrt(x__.error*x__.error+y__.error*y__.error))

  def __rsub__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    return FloatWithError(y__.mean-x__.mean,sqrt(y__.error*y__.error+x__.error*x__.error))

  def __iadd__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    return FloatWithError(x__.mean+y__.mean,sqrt(x__.error*x__.error+y__.error*y__.error))
  def __isub__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    return FloatWithError(x__.mean-y__.mean,sqrt(x__.error*x__.error+y__.error*y__.error))


  def __mul__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    z__mean_      = x__.mean * y__.mean
    #x__rel_error_ = x__.error / x__.mean
    #y__rel_error_ = y__.error / y__.mean
    #return FloatWithError(z__mean_,z__mean_*sqrt(x__rel_error_*x__rel_error_ + y__rel_error_*y__rel_error_))
    z__deri__x__   = y__.mean
    z__deri__y__   = x__.mean
    return FloatWithError(z__mean_,sqrt(z__deri__x__*z__deri__x__*x__.error*x__.error + z__deri__y__*z__deri__y__*y__.error*y__.error))

  def __div__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    z__mean_      = x__.mean / y__.mean
    #x__rel_error_ = x__.error / x__.mean
    #y__rel_error_ = y__.error / y__.mean
    #return FloatWithError(z__mean_,z__mean_*sqrt(x__rel_error_*x__rel_error_ + y__rel_error_*y__rel_error_))
    z__deri__x__   = 1./y__.mean
    z__deri__y__   = -x__.mean/(y__.mean*y__.mean)
    return FloatWithError(z__mean_,sqrt(z__deri__x__*z__deri__x__*x__.error*x__.error + z__deri__y__*z__deri__y__*y__.error*y__.error))

  def __rmul__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    z__mean_      = y__.mean * x__.mean
    #x__rel_error_ = x__.error / x__.mean
    #y__rel_error_ = y__.error / y__.mean
    #return FloatWithError(z__mean_,z__mean_*sqrt(y__rel_error_*y__rel_error_ + x__rel_error_*x__rel_error_))
    z__deri__x__   = y__.mean
    z__deri__y__   = x__.mean
    return FloatWithError(z__mean_,sqrt(z__deri__x__*z__deri__x__*x__.error*x__.error + z__deri__y__*z__deri__y__*y__.error*y__.error))

  def __rdiv__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    z__mean_      = y__.mean / x__.mean
    #x__rel_error_ = x__.error / x__.mean
    #y__rel_error_ = y__.error / y__.mean
    #return FloatWithError(z__mean_,z__mean_*sqrt(y__rel_error_*y__rel_error_ + x__rel_error_*x__rel_error_))
    z__deri__y__   = 1./x__.mean
    z__deri__x__   = -y__.mean/(x__.mean*x__.mean)
    return FloatWithError(z__mean_,sqrt(z__deri__x__*z__deri__x__*x__.error*x__.error + z__deri__y__*z__deri__y__*y__.error*y__.error))

  def __imul__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    z__mean_      = x__.mean * y__.mean
    #x__rel_error_ = x__.error / x__.mean
    #y__rel_error_ = y__.error / y__.mean
    #return FloatWithError(z__mean_,z__mean_*sqrt(x__rel_error_*x__rel_error_ + y__rel_error_*y__rel_error_))
    z__deri__x__   = y__.mean
    z__deri__y__   = x__.mean
    return FloatWithError(z__mean_,sqrt(z__deri__x__*z__deri__x__*x__.error*x__.error + z__deri__y__*z__deri__y__*y__.error*y__.error))

  def __idiv__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    z__mean_      = x__.mean / y__.mean
    #x__rel_error_ = x__.error / x__.mean
    #y__rel_error_ = y__.error / y__.mean
    #return FloatWithError(z__mean_,z__mean_*sqrt(x__rel_error_*x__rel_error_ + y__rel_error_*y__rel_error_))
    z__deri__x__   = 1./y__.mean
    z__deri__y__   = -x__.mean/(y__.mean*y__.mean)
    return FloatWithError(z__mean_,sqrt(z__deri__x__*z__deri__x__*x__.error*x__.error + z__deri__y__*z__deri__y__*y__.error*y__.error))


  def __neg__(self):
    return FloatWithError(0-self.mean,self.error)

  def __pos__(self):
    return FloatWithError(self.mean,self.error)

  def __abs__(self):
    return FloatWithError(abs(self.mean),self.error)


  def __pow__(x__,y__):
    z__mean_      = pow(x__.mean,y__)
    #x__rel_error_ = x__.error / x__.mean
    #return FloatWithError(z__mean_,z__mean_*y__*x__rel_error_)
    z__deri__      = y__*pow(x__.mean,y__-1.)
    return FloatWithError(z__mean_,abs(z__deri__*x__.error))

  def __ipow__(x__,y__):
    z__mean_      = pow(x__.mean,y__)
    #x__rel_error_ = x__.error / x__.mean
    #return FloatWithError(z__mean_,z__mean_*y__*x__rel_error_)
    z__deri__      = y__*pow(x__.mean,y__-1.)
    return FloatWithError(z__mean_,abs(z__deri__*x__.error))


  def sq(self):
    return self**2
  def cb(self):
    return self**3
  def sqrt(self):
    return self**0.5
  def cbrt(self):
    return self**(1./3)


  def exp(self,a__=e):
    y__mean_ = a__**self.mean
    return FloatWithError(y__mean_,abs(y__mean_*log(a__)*self.error))

  def log(self,a__=e):
    y__mean_ = log(self.mean,a__)
    y__deri__ = 1./(self.mean * log(a__))
    return FloatWithError(y__mean_,abs(y__deri__ * self.error))

  
  def sin(self):
    y__mean_ = sin(self.mean)
    y__deri__ = cos(self.mean)
    return FloatWithError(y__mean_,abs(y__deri__ * self.error))

  def cos(self):
    y__mean_ = cos(self.mean)
    y__deri__ = -sin(self.mean)
    return FloatWithError(y__mean_,abs(y__deri__ * self.error))

  def tan(self):
    y__mean_ = tan(self.mean)
    y__deri__ = 1./(cos(self.mean)*cos(self.mean))
    return FloatWithError(y__mean_,abs(y__deri__ * self.error))

  def asin(self):
    y__mean_ = asin(self.mean)
    y__deri__ = 1./sqrt(1. - self.mean*self.mean)
    return FloatWithError(y__mean_,abs(y__deri__ * self.error))

  def acos(self):
    y__mean_ = acos(self.mean)
    y__deri__ = -1./sqrt(1. - self.mean*self.mean)
    return FloatWithError(y__mean_,abs(y__deri__ * self.error))

  def atan(self):
    y__mean_ = atan(self.mean)
    y__deri__ = 1./(1. + self.mean*self.mean)
    return FloatWithError(y__mean_,abs(y__deri__ * self.error))

  
  def sinh(self):
    y__mean_ = sinh(self.mean)
    y__deri__ = cosh(self.mean)
    return FloatWithError(y__mean_,abs(y__deri__ * self.error))

  def cosh(self):
    y__mean_ = cosh(self.mean)
    y__deri__ = sinh(self.mean)
    return FloatWithError(y__mean_,abs(y__deri__ * self.error))

  def tanh(self):
    y__mean_ = tanh(self.mean)
    y__deri__ = 1./(cosh(self.mean)*cosh(self.mean))
    return FloatWithError(y__mean_,abs(y__deri__ * self.error))

  def asinh(self):
    y__mean_ = asinh(self.mean)
    y__deri__ = 1./sqrt(self.mean*self.mean + 1.)
    return FloatWithError(y__mean_,abs(y__deri__ * self.error))

  def acosh(self):
    y__mean_ = acosh(self.mean)
    y__deri__ = 1./sqrt(self.mean*self.mean - 1.)
    return FloatWithError(y__mean_,abs(y__deri__ * self.error))

  def atanh(self):
    y__mean_ = atanh(self.mean)
    y__deri__ = 1./(1. - self.mean*self.mean)
    return FloatWithError(y__mean_,abs(y__deri__ * self.error))
 
