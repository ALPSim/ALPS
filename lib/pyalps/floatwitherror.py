from __future__ import absolute_import
# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009 by Ping Nang (Tama) Ma <pingnang@phys.ethz.ch>
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
from sys  import stdin
from .math import *
        
def get_mean(f):
    try:
        try:
            return f.mean()
        except Exception:
            pass
        return f.mean
    except Exception:
        return f

def get_error(f):
    try:
        return f.error
    except AttributeError:
        return 0

class FloatWithError:


  def __init__(self,mean_=0,error_=0,jackknife=[],binsize=0,timeseries=[]):
    self.mean  = mean_
    self.error = error_
    self.jackknife = jackknife
    self.binsize = binsize
    self.timeseries = list(timeseries)
    try:
      self.shape = self.mean.shape
    except AttributeError as ValueError:
      pass

  def __str__(self):
    return str(self.mean) + ' +/- ' + str(self.error)
  def __expr__(self):
    return expr(self.mean) + ' +/- ' + expr(self.error)
  def __repr__(self):
    return self.__str__()
  
  def __len__(self):
    return len(self.mean)
  
  def __getitem__(self,key):
    try:    jkey = (slice(None),)+key
    except: jkey = (slice(None),key)
    try:    j = self.jackknife[jkey]
    except: j = []
    return FloatWithError(self.mean[key],self.error[key],j)
  
  def __setitem__(self,key,value):
    self.mean [key] = value.mean
    self.error[key] = value.error
    try:    jkey = (slice(None),)+key
    except: jkey = (slice(None),key)
    try:    self.jackknife[jkey] = value.jackknife
    except: self.jackknife = []

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


  def exp(self):
    y__mean_ = exp(self.mean)
    return FloatWithError(y__mean_,abs(y__mean_*self.error))

  def log(self):
    y__mean_ = log(self.mean)
    y__deri__ = 1./(self.mean)
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
 
