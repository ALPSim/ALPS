# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Lukas Gamper <gamperl@gmail.com>
#                       Matthias Troyer <troyer@itp.phys.ethz.ch>
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import pyalps.hdf5 as h5
import numpy as np
import os
import tempfile

def write(ar):
    ar["/int"] =  9
    ar["/double"] =  9.123
    ar["/cplx"] =  complex(1, 2)
    ar["/str"] =  "test"
    ar["/np/int"] =  np.array([1, 2, 3])
    ar["/np2/int"] =  np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    ar["/np/cplx"] =  np.array([[1 + 1j,2 +2j ],[3 + 3j,4 + 4j]])
    
    ar.create_group("/my/group")
    ar["/my/double"] = 9.123
    
    ar.delete_group("/my/group")
    ar.delete_data("/my/double")

def read(ar):
    childs = ar.list_children('/')
    i = ar["/int"]
    d = ar["/double"]
    c = ar["/cplx"]
    s = ar["/str"]
    n = ar["/np/int"]
    x = ar["/np/cplx"]
    
    if len(childs) != 7:
        raise Exception('invalid length of \'/\'')
    if len(ar.extent("/int")) != 1 or ar.extent("/int")[0] != 1 or len(ar.extent("/cplx")) != 1 or ar.extent("/cplx")[0] != 1:
        raise Exception('invalid scalar extent')
    if len(ar.extent("/np/int")) != 1 or ar.extent("/cplx")[0] != 1 or len(ar.extent("/np/cplx")) != 2 or ar.extent("/np/cplx")[0] != 2 or ar.extent("/np/cplx")[1] != 2:
        raise Exception('invalid array extent')
    if type(i) != int or type(d) != float or type(c) != complex or type(s) != str:
        raise Exception('invalid type')
    if i != 9 or d - 9.123 > 0.001 or s != "test" or np.any(n != np.array([1, 2, 3])):
        raise Exception('invalid scalar value')
        
    if np.any(x[0] != np.array([1 + 1j,2 +2j])) or np.any(x[1] != np.array([3 + 3j,4 + 4j])):
        raise Exception('invalid array value')

def test_hdf5():
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "py.h5")
        oar = h5.archive(path, 'w')
        write(oar)
        del oar

        iar = h5.archive(path, 'r')
        if iar.is_complex("/int") or not iar.is_complex("/cplx") or not iar.extent("/np/cplx"):
            raise Exception('invalid complex detection')
        read(iar)
        del iar

        ar = h5.archive(path, 'w')
        write(ar)
        read(ar)
        del ar
    
    print("SUCCESS")


def test_archive_open_state():
    """is_open is a property, and `closed` therefore actually works.

    The legacy Boost.Python module bound the name twice -- add_property
    followed by .def -- and add_to_namespace only merges with an existing
    *function*, so the .def overwrote the property and `ar.is_open` was a bound
    method: always truthy, which left pyalps.hdf5.archive.closed permanently
    False. The nanobind port binds is_open only as a read-only property.
    """
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "open_state.h5")

        ar = h5.archive(path, 'w')
        assert ar.is_open is True
        assert ar.closed is False
        # regression guard: a property, not a callable
        try:
            ar.is_open()
        except TypeError:
            pass
        else:
            raise AssertionError("is_open must be a property, not a method")

        ar['/a'] = 1
        ar.close()
        assert ar.is_open is False
        assert ar.closed is True

        # xml() consults `closed`, so it must now refuse a closed archive
        try:
            ar.xml()
        except h5.ArchiveClosed:
            pass
        else:
            raise AssertionError("xml() on a closed archive must raise ArchiveClosed")

        # the context manager closes on exit
        with h5.archive(path, 'r') as reader:
            assert reader.is_open is True and reader.closed is False
        assert reader.closed is True
