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
