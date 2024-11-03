from __future__ import print_function
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