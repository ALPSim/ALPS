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
  def __init__(self,mean__=0,error__=0):
    self.mean_  = mean__
    self.error_ = error__

  def __str__(self):
    return str(self.mean_) + '\t' + str(self.error_)
  def __repr__(self):
	return self.__str__()
  def __expr__(self):
    return expr(self.mean_) + '\t' + expr(self.error_)
  def mean(self):
    return str(self.mean_)
  def error(self):
    return str(self.error_)

  
  def __add__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    return FloatWithError(x__.mean_+y__.mean_,sqrt(x__.error_*x__.error_+y__.error_*y__.error_))

  def __radd__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    return FloatWithError(y__.mean_+x__.mean_,sqrt(y__.error_*y__.error_+x__.error_*x__.error_))

  def __sub__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    return FloatWithError(x__.mean_-y__.mean_,sqrt(x__.error_*x__.error_+y__.error_*y__.error_))

  def __rsub__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    return FloatWithError(y__.mean_-x__.mean_,sqrt(y__.error_*y__.error_+x__.error_*x__.error_))

  def __iadd__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    return FloatWithError(x__.mean_+y__.mean_,sqrt(x__.error_*x__.error_+y__.error_*y__.error_))
  def __isub__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    return FloatWithError(x__.mean_-y__.mean_,sqrt(x__.error_*x__.error_+y__.error_*y__.error_))


  def __mul__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    z__mean__      = x__.mean_ * y__.mean_
    #x__rel_error__ = x__.error_ / x__.mean_
    #y__rel_error__ = y__.error_ / y__.mean_
    #return FloatWithError(z__mean__,z__mean__*sqrt(x__rel_error__*x__rel_error__ + y__rel_error__*y__rel_error__))
    z__deri__x__   = y__.mean_
    z__deri__y__   = x__.mean_
    return FloatWithError(z__mean__,sqrt(z__deri__x__*z__deri__x__*x__.error_*x__.error_ + z__deri__y__*z__deri__y__*y__.error_*y__.error_))

  def __div__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    z__mean__      = x__.mean_ / y__.mean_
    #x__rel_error__ = x__.error_ / x__.mean_
    #y__rel_error__ = y__.error_ / y__.mean_
    #return FloatWithError(z__mean__,z__mean__*sqrt(x__rel_error__*x__rel_error__ + y__rel_error__*y__rel_error__))
    z__deri__x__   = 1./y__.mean_
    z__deri__y__   = -x__.mean_/(y__.mean_*y__.mean_)
    return FloatWithError(z__mean__,sqrt(z__deri__x__*z__deri__x__*x__.error_*x__.error_ + z__deri__y__*z__deri__y__*y__.error_*y__.error_))

  def __rmul__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    z__mean__      = y__.mean_ * x__.mean_
    #x__rel_error__ = x__.error_ / x__.mean_
    #y__rel_error__ = y__.error_ / y__.mean_
    #return FloatWithError(z__mean__,z__mean__*sqrt(y__rel_error__*y__rel_error__ + x__rel_error__*x__rel_error__))
    z__deri__x__   = y__.mean_
    z__deri__y__   = x__.mean_
    return FloatWithError(z__mean__,sqrt(z__deri__x__*z__deri__x__*x__.error_*x__.error_ + z__deri__y__*z__deri__y__*y__.error_*y__.error_))

  def __rdiv__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    z__mean__      = y__.mean_ / x__.mean_
    #x__rel_error__ = x__.error_ / x__.mean_
    #y__rel_error__ = y__.error_ / y__.mean_
    #return FloatWithError(z__mean__,z__mean__*sqrt(y__rel_error__*y__rel_error__ + x__rel_error__*x__rel_error__))
    z__deri__y__   = 1./x__.mean_
    z__deri__x__   = -y__.mean_/(x__.mean_*x__.mean_)
    return FloatWithError(z__mean__,sqrt(z__deri__x__*z__deri__x__*x__.error_*x__.error_ + z__deri__y__*z__deri__y__*y__.error_*y__.error_))

  def __imul__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    z__mean__      = x__.mean_ * y__.mean_
    #x__rel_error__ = x__.error_ / x__.mean_
    #y__rel_error__ = y__.error_ / y__.mean_
    #return FloatWithError(z__mean__,z__mean__*sqrt(x__rel_error__*x__rel_error__ + y__rel_error__*y__rel_error__))
    z__deri__x__   = y__.mean_
    z__deri__y__   = x__.mean_
    return FloatWithError(z__mean__,sqrt(z__deri__x__*z__deri__x__*x__.error_*x__.error_ + z__deri__y__*z__deri__y__*y__.error_*y__.error_))

  def __idiv__(x__,y__):
    if (isinstance(y__,float) | isinstance(y__,int)):
      y__ = FloatWithError(y__)
    z__mean__      = x__.mean_ / y__.mean_
    #x__rel_error__ = x__.error_ / x__.mean_
    #y__rel_error__ = y__.error_ / y__.mean_
    #return FloatWithError(z__mean__,z__mean__*sqrt(x__rel_error__*x__rel_error__ + y__rel_error__*y__rel_error__))
    z__deri__x__   = 1./y__.mean_
    z__deri__y__   = -x__.mean_/(y__.mean_*y__.mean_)
    return FloatWithError(z__mean__,sqrt(z__deri__x__*z__deri__x__*x__.error_*x__.error_ + z__deri__y__*z__deri__y__*y__.error_*y__.error_))


  def __neg__(self):
    return FloatWithError(0-self.mean_,self.error_)

  def __pos__(self):
    return FloatWithError(self.mean_,self.error_)

  def __abs__(self):
    return FloatWithError(abs(self.mean_),self.error_)


  def __pow__(x__,y__):
    z__mean__      = pow(x__.mean_,y__)
    #x__rel_error__ = x__.error_ / x__.mean_
    #return FloatWithError(z__mean__,z__mean__*y__*x__rel_error__)
    z__deri__      = y__*pow(x__.mean_,y__-1.)
    return FloatWithError(z__mean__,abs(z__deri__*x__.error_))

  def __ipow__(x__,y__):
    z__mean__      = pow(x__.mean_,y__)
    #x__rel_error__ = x__.error_ / x__.mean_
    #return FloatWithError(z__mean__,z__mean__*y__*x__rel_error__)
    z__deri__      = y__*pow(x__.mean_,y__-1.)
    return FloatWithError(z__mean__,abs(z__deri__*x__.error_))


  def sq(self):
    return self**2
  def cb(self):
    return self**3
  def sqrt(self):
    return self**0.5
  def cbrt(self):
    return self**(1./3)


  def exp(self,a__=e):
    y__mean__ = a__**self.mean_
    return FloatWithError(y__mean__,abs(y__mean__*log(a__)*self.error_))

  def log(self,a__=e):
    y__mean__ = log(self.mean_,a__)
    y__deri__ = 1./(self.mean_ * log(a__))
    return FloatWithError(y__mean__,abs(y__deri__ * self.error_))

  
  def sin(self):
    y__mean__ = sin(self.mean_)
    y__deri__ = cos(self.mean_)
    return FloatWithError(y__mean__,abs(y__deri__ * self.error_))

  def cos(self):
    y__mean__ = cos(self.mean_)
    y__deri__ = -sin(self.mean_)
    return FloatWithError(y__mean__,abs(y__deri__ * self.error_))

  def tan(self):
    y__mean__ = tan(self.mean_)
    y__deri__ = 1./(cos(self.mean_)*cos(self.mean_))
    return FloatWithError(y__mean__,abs(y__deri__ * self.error_))

  def asin(self):
    y__mean__ = asin(self.mean_)
    y__deri__ = 1./sqrt(1. - self.mean_*self.mean_)
    return FloatWithError(y__mean__,abs(y__deri__ * self.error_))

  def acos(self):
    y__mean__ = acos(self.mean_)
    y__deri__ = -1./sqrt(1. - self.mean_*self.mean_)
    return FloatWithError(y__mean__,abs(y__deri__ * self.error_))

  def atan(self):
    y__mean__ = atan(self.mean_)
    y__deri__ = 1./(1. + self.mean_*self.mean_)
    return FloatWithError(y__mean__,abs(y__deri__ * self.error_))

  
  def sinh(self):
    y__mean__ = sinh(self.mean_)
    y__deri__ = cosh(self.mean_)
    return FloatWithError(y__mean__,abs(y__deri__ * self.error_))

  def cosh(self):
    y__mean__ = cosh(self.mean_)
    y__deri__ = sinh(self.mean_)
    return FloatWithError(y__mean__,abs(y__deri__ * self.error_))

  def tanh(self):
    y__mean__ = tanh(self.mean_)
    y__deri__ = 1./(cosh(self.mean_)*cosh(self.mean_))
    return FloatWithError(y__mean__,abs(y__deri__ * self.error_))

  def asinh(self):
    y__mean__ = asinh(self.mean_)
    y__deri__ = 1./sqrt(self.mean_*self.mean_ + 1.)
    return FloatWithError(y__mean__,abs(y__deri__ * self.error_))

  def acosh(self):
    y__mean__ = acosh(self.mean_)
    y__deri__ = 1./sqrt(self.mean_*self.mean_ - 1.)
    return FloatWithError(y__mean__,abs(y__deri__ * self.error_))

  def atanh(self):
    y__mean__ = atanh(self.mean_)
    y__deri__ = 1./(1. - self.mean_*self.mean_)
    return FloatWithError(y__mean__,abs(y__deri__ * self.error_))
 
