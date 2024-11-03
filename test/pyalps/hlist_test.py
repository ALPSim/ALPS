# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Ping Nang Ma
#
# This software is part of the ALPS libraries, published under the ALPS
# Library License; you can use, redistribute it and/or modify it under
# the terms of the license, either version 1 or (at your option) any later
# version.
#  
# You should have received a copy of the ALPS Library License along with
# the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
# available from http://alps.comp-phys.org/.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
# FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
# 
# ****************************************************************************
import pytest

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
