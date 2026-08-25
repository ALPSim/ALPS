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

import numpy as np

def _values_equal(left, right):
    try:
        return bool(np.all(left == right))
    except (TypeError, ValueError):
        return False

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
            if not _values_equal(val0, idict[key]):
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
            if not _values_equal(val0, idict[key]):
                take = False
        if not take:
            ret.append(key)
    return ret
