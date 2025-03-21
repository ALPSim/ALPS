# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Olivier Parcollet
# 
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
# 
# ****************************************************************************

""" 
This module contains the archive class, used to read / write hdf5 files
"""

from .cxx.pyngshdf5_c import hdf5_archive_impl
from .cxx.pyngshdf5_c import register_archive_exception_type

class ArchiveError(Exception): pass
register_archive_exception_type(0, ArchiveError)

class ArchiveNotFound(ArchiveError, IOError): pass
register_archive_exception_type(1, ArchiveNotFound)

class ArchiveClosed(ArchiveError, ValueError): pass
register_archive_exception_type(2, ArchiveClosed)

class InvalidPath(ArchiveError, SyntaxError): pass
register_archive_exception_type(3, InvalidPath)

class PathNotFound(ArchiveError, LookupError): pass
register_archive_exception_type(4, PathNotFound)

class WrongType(ArchiveError, TypeError): pass
register_archive_exception_type(5, WrongType)

del register_archive_exception_type

#TODO: move to hdf5 module
class archive(hdf5_archive_impl):
    def __init__(self, filename, mode = 'r'):
        hdf5_archive_impl.__init__(self, filename, mode)

    @property # TODO: move to c++
    def closed(self):
        return not self.is_open

    def __enter__(self, *args, **kwargs): # TODO: move to c++
        return self

    def __exit__(self, *args, **kwargs): # TODO: move to c++
        self.close()

    def xml(self, path="/"):
        """
        Returns an XML formatted string of the archive.
        This function still needs to have attributes implemented.
        """
        if self.closed: raise ArchiveClosed("I/O operation on closed file")
        if self.is_group(path):
            ret = ""
            for child in self.list_children(path):
                ret += "<%s>\n" % child
                ret += self.as_xml(path+("" if path[-1] == "/" else "/" )+child)
                ret += "\n</%s>" % child
            return ret
        else:
            return "%s" %(str(self[path]))

# TODO: remove this in future
import warnings

def h5ar(*args, **kwargs):
    warnings.warn("The object 'h5ar' is deprecated and will be removed. Use 'archive' instead.", DeprecationWarning)
    return archive(*args, **kwargs)

class iArchive(archive):
    def __init__(self, filename):
        warnings.warn("The object 'iArchive' is deprecated and will be removed. Use 'archive' instead.", DeprecationWarning)
        archive.__init__(self, str(filename), 'r')
    def read(self, path):
        return archive.__getitem__(self, path)

class oArchive(archive):
    def __init__(self, filename):
        warnings.warn("The object 'oArchive' is deprecated and will be removed. Use 'archive' instead.", DeprecationWarning)
        archive.__init__(self, str(filename), 'w')
    def read(self, path):
        return archive.__getitem__(self, path)
    def write(self, path, value):
        archive.__setitem__(self, path, value)
