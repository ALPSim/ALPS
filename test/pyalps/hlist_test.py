# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Ping Nang Ma
#
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************
from pyalps.hlist import HList


def test_hlist():
    
    hl = HList([[1,2,3],[4,5]])
    
    print(hl)
    assert list(hl) == [1,2,3,4,5]
    # [[1,2,3],[4,5]]
    
    # !!! Testing linear access
    
    print(hl[0])
    assert hl[0] == 1
    # 1
    
    print(hl[0:2])
    assert hl[0:2] == [1, 2]
    # [1, 2]
    
    # !!! Testing 'recursive' access
    
    print(hl[0,0])
    assert hl[0,0] == 1
    # 1
    print(hl[1,1])
    assert hl[1,1] == 5
    # 5
    
    # !!! Linear assignment
    hl[0] = 27
    print(hl[0])
    assert hl[0] == 27
    print(hl[0,0])
    assert hl[0,0] == 27
    # 27
    
    hl[1,1] = 13
    print(hl[1,1])
    assert hl[1,1] == 13 
    print(hl[4])
    assert hl[4] == 13
    # 13


if __name__ == '__main__':
    test_hlist()
