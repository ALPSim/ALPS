 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 #                                                                                 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   #
 #                                                                                 #
 # ALPS Project: https://alps.comp-phys.org/                                       #
 # SPDX-License-Identifier: MIT                                                    #
 #                                                                                 #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import pyalps.hdf5 as hdf5
import pyalps.ngs as ngs
import os
import tempfile

orig_dict = {
    'val1' : 42,
    'val2' : '42',
    'a' : 1,
    'x' : 2,
    'b' : 3
}
# Values are supplied as text, as an INI file supplies them. Until a program
# defines a parameter's type, it reads back as a number if its text is one,
# so the string '42' reads back as the integer 42.
expected = dict(orig_dict, val2=42)
def assert_type(p, k):
    assert type(p[k]) == type(expected[k])

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
        assert p[k] == expected[k]
        assert p.supplied(k)
        assert_type(p, k)
        print(k,'ok!')
    ## Check nonetype
    assert type(p["undefined"]) == type(None)
    
    ## Write to and load from hdf5 without leaving test artifacts in the tree.
    with tempfile.TemporaryDirectory() as directory:
        parms1 = os.path.join(directory, 'parms1.h5')
        parms2 = os.path.join(directory, 'parms2.h5')
        with hdf5.archive(parms1, 'w') as oar:
            p.save(oar) # does not use path '/parameters'

        with hdf5.archive(parms2, 'w') as oar:
            for key in sorted(p.keys()):
                print(key)
                oar['parameters/' + key] = p[key]

        # Preserve the existing simultaneous-reader exercise.
        with hdf5.archive(parms2, 'r'):
            with hdf5.archive(parms2, 'r') as iar:
                p.load(iar)

        for k in sorted(orig_dict.keys()):
            assert p[k] == expected[k]
            assert_type(p, k)
            print(k,'ok!')

if __name__ == '__main__':
    test_params()
