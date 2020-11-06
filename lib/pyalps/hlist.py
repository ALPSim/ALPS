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

def depth(hl):
    ret = 0
    if len(hl) == 0:
        return 0
    elif type(hl[0]) == list:
        return 1+depth(hl[0])
    else:
        return 1

def index__(hl,indices,idx,level,max_level):
    for ie in range(len(hl)):
        idx[level] = ie
        e = hl[ie]
        if level < max_level-1:
            index__(e,indices,idx,level+1,max_level)
        else:
            indices.append(tuple(idx))

def hset__(hl,idx,level,value):
    if level == len(idx)-1:
        hl[idx[level]] = value
    else:
        hset__(hl[idx[level]],idx,level+1,value)

def hget__(hl,idx,level):
    if level == len(idx)-1:
        return hl[idx[level]]
    else:
        return hget__(hl[idx[level]],idx,level+1)

def flatten(sl, fdepth = None):
    """ turns a hierarchical list of lists into a flat list 
    
        this function turns a hierarchical list (a list of lists of lists ....) into just a flat list.
        The parameters are the hierachical list, and optionally a depth at which the flattening should happen, to keep, e.g. the top-most structure.
    """
    
    if fdepth == None:
        fdepth = depth(list(sl))
    return HList(list(sl), fdepth)

def deep_flatten(sl, fdepth = None):
    return [x for x in flatten(sl, fdepth)]

def copy_structure(sl):
    try:
        ll = len(sl)
        return [copy_structure(isl) for isl in sl]
    except Exception:
        return 1

def happly(functor, sl, fdepth = None, params = None):
    if fdepth == None:
        fdepth = depth(sl)
    hl = HList(sl, fdepth)
    hl.apply(functor, params)

def hmap(functor, sl, fdepth = None, params = None):
    if fdepth == None:
        fdepth = depth(sl)
    hl = HList(sl, fdepth)
    return hl.map(functor, params)

class HList:
    def __init__(self):
        self.data_ = []
        self.indices_ = []
    
    def __init__(self,init,fdepth = None):
        self.data_ = init
        self.indices_ = []
        
        if fdepth == None:
            fdepth = depth(self.data_)
        if fdepth < 0:
            fdepth = depth(self.data_) + fdepth
        
        index__(self.data_, self.indices_, [0 for q in range(depth(self.data_))], 0, fdepth)
        self.indices_ = [idx[0:fdepth] for idx in self.indices_]
    
    def __len__(self):
        return len(self.indices_)
    
    def __getitem__(self, key):
        if type(key) == tuple:
            return hget__(self.data_,key,0)
        elif type(key) == list:
            return [self[k] for k in key]
        else:
            return self[self.indices_[key]]
    
    def __repr__(self):
        return str([self[k] for k in self.indices_])
    
    def __setitem__(self, key, value):
        if type(key) == tuple:
            hset__(self.data_,key,0,value)
        elif type(key) == list:
            raise TypeError("Assigning to slices is not supported")
        else:
            self[self.indices_[key]] = value
    
    def indices(self):
        return self.indices_
    
    def data(self):
        return self.data_
    
    def apply(self, functor, params = None):
        for idx in self.indices_:
            if params == None:
                self[idx] = functor(self[idx])
            else:
                self[idx] = functor(self[idx], params)
    
    def map(self, functor, params = None):
        ret = copy_structure(self.data_)
        rethl = HList(ret)
        for idx in self.indices_:
            if params == None:
                rethl[idx] = functor(self[idx])
            else:
                rethl[idx] = functor(self[idx], params)
        return ret

def hlist_to_dict(hl, key = 'Observable'):
    f = lambda data: dict( [(q.props[key], q) for q in data] )
    return hmap(f, hl, -1)

