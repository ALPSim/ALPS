# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 1994-2009 by Bela Bauer <bauerb@phys.ethz.ch>
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


