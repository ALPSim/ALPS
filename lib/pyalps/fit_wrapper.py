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
from scipy import optimize


__all__ = ['fit', 'Parameter']

class Parameter:    
    def __init__(self, value):
        self.value = value

    def set(self, value):
        self.value = value

    def get(self):
        return self.value

    def __call__(self):
        return self.value

def fit(self,function, parameters, y, x = np.arange(0)):
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return y - function(self,x,parameters)

    if len(x) != y.shape[0]:
        x = np.arange(y.shape[0])
    p = [param() for param in parameters]
    optimize.leastsq(f, p)
