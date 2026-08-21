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
# methods onto the concrete classes instead. A method is copied when the
# class doesn't provide its own — "inherited from object" counts as absent,
# otherwise __eq__/__ne__ (which every type inherits from object) would be
# skipped and mapping equality lost. __hash__ rides along as None, exactly
# as MutableMapping inheritance made these types unhashable before.
for _mapping_type in (params, observables, results):
    MutableMapping.register(_mapping_type)
    for _method in ("keys", "values", "items", "get", "pop", "popitem",
                    "clear", "update", "setdefault", "__eq__", "__ne__",
                    "__hash__"):
        if getattr(_mapping_type, _method, None) is getattr(object, _method, None):
            setattr(_mapping_type, _method, getattr(MutableMapping, _method))

# Two mixin methods cannot simply be copied onto these types.
#
# MutableMapping.pop reads `self.__marker`, which name-mangles to
# self._MutableMapping__marker -- a class attribute of MutableMapping. Under
# the old `params.__bases__ = (MutableMapping,) + ...` rebasing that resolved
# through the MRO; on a virtual subclass whose methods were copied it does
# not, so pop() raised AttributeError instead of returning the default or
# raising KeyError. Supply an implementation that tests membership instead of
# relying on __getitem__ raising.
#
# alps::params compounds this: its __getitem__ returns None for an undefined
# key rather than raising KeyError (inherited from the Boost.Python module and
# pinned by test/pyalps/pyparams_test.py), so get() ignored its default and
# setdefault() returned None while storing nothing.
_MAPPING_POP_MARKER = object()


def _mapping_pop(self, key, default=_MAPPING_POP_MARKER):
    if key in self:
        value = self[key]
        del self[key]
        return value
    if default is _MAPPING_POP_MARKER:
        raise KeyError(key)
    return default


def _params_get(self, key, default=None):
    return self[key] if key in self else default


def _params_setdefault(self, key, default=None):
    if key in self:
        return self[key]
    self[key] = default
    return default


for _mapping_type in (params, observables, results):
    _mapping_type.pop = _mapping_pop

params.get = _params_get
params.setdefault = _params_setdefault

from .cxx.pyngsbase_c import mcbase

from .cxx.pyngsapi_c import collectResults, saveResults

from .cxx.pyngsrandom01_c import random01
