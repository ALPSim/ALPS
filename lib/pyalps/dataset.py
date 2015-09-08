# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 1994-2009 by Bela Bauer <bauerb@phys.ethz.ch>
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

import copy
import numpy as np

class ResultProperties:
    def __init__(self):
        self.props = {}

class DataSet(ResultProperties):
    """
    The DataSet class stores a set of data, usually in XY format, along with all the properties
    describing the data, such as input parameters to the simulation etc.
    
    Members are:
     * x, y - These contain the data and are expected to come as lists of Numpy arrays
              by many functions operating on DataSets. However, for user-supplied functions,
              other ways of representing data may be used.
     * props - This is a dictionary of properties describing the dataset.
    """
    def __init__(self,x=None,y=None,props=None):
        ResultProperties.__init__(self)
        if x is None:   self.x = np.array([])
        else:           self.x = x
        if y is None:   self.y = np.array([])
        else:           self.y = y
        if props is not None:   self.props = props
    
    def __repr__(self):
        return "x=%s\ny=%s\nprops=%s" % (self.x, self.y, self.props)
        
class ResultFile(ResultProperties):
    def __init__(self,fn=None):
        ResultProperties.__init__(self)
        if fn is not None:
            self.props['filename'] = fn


