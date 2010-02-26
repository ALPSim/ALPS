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
    def __deepcopy__(self,memo):
        ret = ResultProperties()
        ret.props = copy.deepcopy(self.props,memo)
        return ret
    
    def __init__(self):
        self.props = {}

class DataSet(ResultProperties):
    def __deepcopy__(self,memo):
        ret = DataSet()
        ret.props = copy.deepcopy(self.props,memo)
        ret.x = copy.deepcopy(self.x,memo)
        ret.y = copy.deepcopy(self.y,memo)
        return ret

    def __init__(self):
        ResultProperties.__init__(self)
        
        self.x = np.array([])
        self.y = np.array([])
        
class ResultFile(ResultProperties):
    def __deepcopy__(self,memo):
        ret = ResultFile
        ret.props = copy.deepcopy(self.props,memo)
        return ret
    
    def __init__(self,fn=None):
        ResultProperties.__init__(self)
        if fn != None:
            self.props['filename'] = fn
        