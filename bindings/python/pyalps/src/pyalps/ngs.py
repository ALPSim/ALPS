 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 #                                                                                 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   #
 #                      2012 by Troels F. Roennow <tfr@nanophysics.dk>             #
 #                                                                                 #
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
 #                                                                                 #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from collections.abc import MutableMapping
from .cxx.pyngsparams_c import params

from .cxx.pyngsobservable_c import observable
def _observable_lshift(self, other):
    self.append(other)
    return self
observable.__lshift__ = _observable_lshift

class RealObservable:
    def __init__(self, name, binnum = 0):
        self.name = name
        self.binnum = binnum
    def addToObservables(self, observables): #rename this with new ALEA
        observables.createRealObservable(self.name, self.binnum)

class RealVectorObservable:
    def __init__(self, name, binnum = 0):
        self.name = name
        self.binnum = binnum
    def addToObservables(self, observables): #rename this with new ALEA
        observables.createRealVectorObservable(self.name, self.binnum)

from .cxx.pyngsobservables_c import observables

from .cxx.pyngsobservable_c import createRealObservable #remove this with new ALEA!
from .cxx.pyngsobservable_c import createRealVectorObservable #remove this with new ALEA!

from .cxx.pyngsresult_c import result
from .cxx.pyngsresult_c import observable2result #remove this with new ALEA!

from .cxx.pyngsresults_c import results

# Boost.Python allowed mutating extension-type base classes after creation.
# nanobind extension types use a different allocator/deallocator layout, so
# register them as virtual MutableMapping implementations and copy the mixin
# methods onto the concrete classes instead.
for _mapping_type in (params, observables, results):
    MutableMapping.register(_mapping_type)
    for _method in ("keys", "values", "items", "get", "pop", "popitem",
                    "clear", "update", "setdefault", "__eq__", "__ne__"):
        if not hasattr(_mapping_type, _method):
            setattr(_mapping_type, _method, getattr(MutableMapping, _method))

from .cxx.pyngsbase_c import mcbase

from .cxx.pyngsapi_c import collectResults, saveResults

from .cxx.pyngsrandom01_c import random01
