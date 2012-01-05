 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 #                                                                                 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   #
 #                                                                                 #
 # This software is part of the ALPS libraries, published under the ALPS           #
 # Library License; you can use, redistribute it and/or modify it under            #
 # the terms of the license, either version 1 or (at your option) any later        #
 # version.                                                                        #
 #                                                                                 #
 # You should have received a copy of the ALPS Library License along with          #
 # the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       #
 # available from http://alps.comp-phys.org/.                                      #
 #                                                                                 #
 #  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     #
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        #
 # FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       #
 # SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       #
 # FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     #
 # ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     #
 # DEALINGS IN THE SOFTWARE.                                                       #
 #                                                                                 #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from collections import MutableMapping
import types

from pyngsparams_c import *
params.__bases__ = (MutableMapping, ) + params.__bases__

from pyngsobservable_c import *
class ObservableOperators:
    def __lshift__(self, other):
        self.append(other)
observable.__bases__ = (ObservableOperators, ) + observable.__bases__

class RealObservable:
    def __init__(self, name, binnum = 0):
        self.name = name
        self.binnum = binnum
    def addToObservables(self, observables):
        observables.createRealObservable(self.name, self.binnum)

class RealVectorObservable:
    def __init__(self, name, binnum = 0):
        self.name = name
        self.binnum = binnum
    def addToObservables(self, observables):
        observables.createRealVectorObservable(self.name, self.binnum)

from pyngsobservables_c import *
class ObservablesOperators:
    def __lshift__(self, other):
        other.addToObservables(self)
observables.__bases__ = (ObservablesOperators, MutableMapping, ) + observables.__bases__

from pyngsresult_c import *

from pyngsresults_c import *
results.__bases__ = (MutableMapping, ) + results.__bases__

from pyngsbase_c import *
class base(base_impl):
    def run(self, callback = lambda: True):
        base_impl.run(self, callback)

from pyngshdf5_c import *
class h5ar(hdf5_archive_impl):
    READ = 0
    WRITE = 1
    REPLACE = 2
    COMPRESS = 4
    def __init__(self, filename, mode = WRITE):
        hdf5_archive_impl.__init__(self, filename, mode)
    def __getitem__(self, path):
        return self.load(path)
    def __setitem__(self, path, value):
        self.save(path, value)
    def save(self, path, data):
        if hasattr(data, 'save') and type(getattr(data, 'save')) == types.MethodType:
            current = self.context
            self.set_context(path)
            data.save(self)
            self.set_context(current)
        else:
            dtp = type(data).__name__
            if dtp == 'ndarray':
                dtp = str(data.dtype)
            hdf5_archive_impl.save(self, path, data, dtp)

from pyngsapi_c import *
