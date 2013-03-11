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

import numpy as np

def dict_intersect(dicts):
    """ computes the intersection of a list of dicts
    
        this function takes a list of dicts as input and returns a dict containing all those key-value pairs that appear with identical values in all dicts 
    """
    sets = [set(q.keys()) for q in dicts]
    intersection = sets[0]
    for iset in sets:
        intersection &= iset
    ret = {}
    for key in intersection:
        take = True
        val0 = dicts[0][key]
        for idict in dicts:
            try:
                if val0 != idict[key]:
                    take = False
            except:
                if np.all(val0 != idict[key]):
                    take = False
        if take:
            ret[key] = dicts[0][key]
    return ret

def dict_difference(dicts):
    sets = [set(q.keys()) for q in dicts]
    intersection = sets[0]
    for iset in sets:
        intersection &= iset
    ret = []
    for key in intersection:
        take = True
        val0 = dicts[0][key]
        for idict in dicts:
            if val0 != idict[key]:
                take = False
        if not take:
            ret.append(key)
    return ret
