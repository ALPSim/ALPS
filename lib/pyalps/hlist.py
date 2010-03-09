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
    if fdepth == None:
        fdepth = depth(sl)
    return HList(sl, fdepth)

def deep_flatten(sl, fdepth = None):
    return [x for x in flatten(sl, fdepth)]

class HList:
    def __init__(self):
        self.data_ = []
        self.indices_ = []
    
    def __init__(self,init,fdepth = None):
        self.data_ = init
        self.indices_ = []
        
        if fdepth == None:
            fdepth = depth(self.data_)
        
        index__(self.data_, self.indices_, [0 for q in range(depth(self.data_))], 0, fdepth)
        self.indices_ = [idx[0:fdepth] for idx in self.indices_]
    
    def __getitem__(self, key):
#        print 'Get key',key
        if type(key) == tuple:
            return hget__(self.data_,key,0)
        elif type(key) == list:
            return [self[k] for k in key]
        else:
            return self[self.indices_[key]]
    
    def __repr__(self):
        return str([self[k] for k in self.indices_])
    
    def __setitem__(self, key, value):
#        print 'Set key',key
        if type(key) == tuple:
            hset__(self.data_,key,0,value)
        elif type(key) == list:
            raise TypeError("Assigning to slices is not supported")
        else:
            self[self.indices_[key]] = value
    
