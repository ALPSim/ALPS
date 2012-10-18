 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 #                                                                                 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   #
 #                      2012 by Troels F. Roennow <tfr@nanophysics.dk>             #
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
class mcbase(mcbase_impl):
    def run(self, callback = lambda: True):
        mcbase_impl.run(self, callback)

from pyngsapi_c import *

from pyngshdf5_c import hdf5_archive_impl

class ArchiveIOException(IOError):
    pass

# TODO: export hdf5archive errors from c++
class ArchivePathException(ArchiveIOException):
    def __init__(self, msg, path):
        super(ArchivePathException,self).__init__(msg)
        self.path = path

# TODO: either do hdf5.archive or renmae it to hdf5Archive
class archive:
    def __init__(self, filename, mode = "r", *args, **kwargs):
        """
        Opens an archive in specified mode. If no mode is specified 'r'
        is used as default.
        """
        self._archive = hdf5_archive_impl(filename, mode, *args, **kwargs)

        def wrapper(name):
            def f(*args, **kwargs):
                if self.closed: raise ArchiveIOException("I/O operation on closed file")
                return getattr(self._archive, name)(*args, **kwargs)
            return f
        
        excp = [] #["__setitem__", "__getitem__"]
        for name in dir(self._archive):
            if name.startswith("_") and name not in excp: continue
            setattr(self, name, wrapper(name) )

    @property
    def closed(self):
        """
        Property to check whether the file is open or closed.
        """
        return self._archive is None

    def __getitem__(self, path):
        """
        The getter of the archive is a work around for the missing C++
        implementation of the exceptions. 
        TODO: 
        """
        try:
            return self._archive[path]
        except:
            raise ArchivePathException("could not read '%s'"%path, path)

    def __setitem__(self, path, value):
        """
        The setter of the archive is a work around for the missing C++
        implementation of the exceptions.  
        TODO: 
        """
        try:
            self._archive[path] = value
        except:
            raise ArchivePathException("no write access to '%s'"%path, path)

    def close(self):
        """
        Closes the archive.
        """
        del self._archive
        self._archive = None

    def __enter__(self, *args, **kwargs):
        return self

    def __exit__(self, *args, **kwargs):
        self.close()
        return self

    def xml(self,path="/"):
        """
        Returns an XML formatted string of the archive.         
        This function still needs to have attributes implemented.
        """
        if self.closed: ArchiveIOException("I/O operation on closed file")
        if self._archive.is_group(path):
            ret = ""
            for child in self._archive.list_children(path):
                ret += "<%s>\n" % child
                ret += self.as_xml(path+("" if path[-1] == "/" else "/" )+child) 
                ret += "\n</%s>" % child
            return ret
        else:
            return "%s" %(str(self._archive[path]))

import warnings
def h5ar(*args, **kwargs):
    warnings.warn("The object 'h5ar' is deprecated and will be removed. Use 'archive' instead.", DeprecationWarning)
    return archive(*args, **kwargs)
        
