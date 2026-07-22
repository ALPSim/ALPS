 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 #                                                                                 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   #
 #                      2012 by Troels F. Roennow <tfr@nanophysics.dk>             #
 #                                                                                 #
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
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
