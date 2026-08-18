# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 1994-2009 by Bela Bauer <bauerb@phys.ethz.ch>
# 
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
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


