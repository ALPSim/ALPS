 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 #                                                                                 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   #
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

#TODO: werden funktionen als camel notation oder mit unterstrich benannt?
from collections import MutableMapping
import types

from pyngsparams_c import *
params.__bases__ = (MutableMapping, ) + params.__bases__

from pyngsobservable_c import *
observable.__bases__ = (MutableMapping, ) + observable.__bases__

from pyngsobservables_c import *
observables.__bases__ = (MutableMapping, ) + observables.__bases__

from pyngsresult_c import *
result.__bases__ = (MutableMapping, ) + result.__bases__

from pyngsresults_c import *
results.__bases__ = (MutableMapping, ) + results.__bases__

from pyngsbase_c import *
class base(base_impl):
    def __init__(self, arg = {}):
        if type(arg).__name__ == "dict":
            self.params = params(arg)
        else:
            self.params = arg
        base_impl.__init__(self, arg)
        self.measurements = self.get_measurements()
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
            while dtp == 'ndarray':
                dtp = type(el = data[0]).__name__
            hdf5_archive_impl.save(self, path, data, dtp)


from pyngsapi_c import *
