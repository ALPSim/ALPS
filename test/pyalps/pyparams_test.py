from __future__ import print_function
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 #                                                                                 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   #
 #                                                                                 #
 # SPDX-License-Identifier: MIT
 #                                                                                 #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import pyalps.hdf5 as hdf5
import pyalps.ngs as ngs
import sys

orig_dict = {
    'val1' : 42,
    'val2' : '42',
    'a' : 1,
    'x' : 2,
    'b' : 3
}
def assert_type(p, k):
    assert type(p[k]) == type(orig_dict[k])

def test_params():    
    ## Create params
    p = ngs.params({
        'val1' : 42,
        'val2' : '42',
        'a' : 1,
        'x' : 2,
        'b' : 3
    })
    ## check content
    for k in sorted(orig_dict.keys()):
        assert p[k] == orig_dict[k]
        assert_type(p, k)
        print(k,'ok!')
    ## Check nonetype
    assert type(p["undefined"]) == type(None)
    
    ## Write to hdf5
    with hdf5.archive('parms1.h5', 'w') as oar:
        p.save(oar) # does not use path '/parameters'
    
    with hdf5.archive('parms2.h5', 'w') as oar:
        for key in sorted(p.keys()):
            print(key)
            oar['parameters/' + key] = p[key]
    ## Load from hdf5
    with hdf5.archive('parms2.h5', 'r') as oar:
        iar = hdf5.archive('parms2.h5', 'r')
        p.load(iar)
    
        for k in sorted(orig_dict.keys()):
            assert p[k] == orig_dict[k]
            assert_type(p, k)
            print(k,'ok!')

if __name__ == '__main__':
    test_params()