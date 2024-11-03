from __future__ import print_function
# ****************************************************************************
#
# ALPS Project: Algorithms and Libraries for Physics Simulations
#
# ALPS Libraries
#
# Copyright (C) 2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch> ,
#                       Lukas Gamper <gamperl@gmail.com>
#                       Matthias Troyer <troyer@itp.phys.ethz.ch>
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

from pyalps.alea import *
import numpy as np

def str_prec(a):
    return '{:.12f}'.format(a)


def test_mcdata():

    print("\nTesting MCScalarData")
    print("\n------------------------\n")
    
    a = MCScalarData(0.81,0.1)
    b = MCScalarData(1.21,0.15)
    c = MCScalarData(-1.5,0.2)
    
    print("Initialization:\n")
    print("a:\t" + str(a))
    print("b:\t" + str(b))
    print("c:\t" + str(c))
    
    print("\n")
    
    print("Operation:\n")
    
    a += b
    print("a += b:\t" + str_prec(a))
    a = MCScalarData(1.2,0.1)
    
    a -= b
    print("a -= b:\t" + str_prec(a))
    a = MCScalarData(1.2,0.1)
    
    a *= b
    print("a *= b:\t" + str_prec(a))
    a = MCScalarData(1.2,0.1)
    
    a /= b
    print("a /= b:\t" + str_prec(a))
    a = MCScalarData(1.2,0.1)
    
    print("\n")
    
    a += 2.
    print("a += 2.:\t" + str_prec(a))
    a = MCScalarData(1.2,0.1)
    
    a -= 2.
    print("a -= 2.:\t" + str_prec(a))
    a = MCScalarData(1.2,0.1)
    
    a *= 2.
    print("a *= 2.:\t" + str_prec(a))
    a = MCScalarData(1.2,0.1)
    
    a /= 2.
    print("a /= 2.:\t" + str_prec(a))
    a = MCScalarData(1.2,0.1)
    
    print("\n")
    
    print("a + b:\t" + str_prec(a / b))
    print("a + 2.:\t" + str_prec(a / 2.))
    print("2. + a:\t" + str_prec(2. / a))
    print("a - b:\t" + str_prec(a / b))
    print("a - 2.:\t" + str_prec(a / 2.))
    print("2. - a:\t" + str_prec(2. / a))
    print("a * b:\t" + str_prec(a / b))
    print("a * 2.:\t" + str_prec(a / 2.))
    print("2. * a:\t" + str_prec(2. / a))
    print("a / b:\t" + str_prec(a / b))
    print("a / 2.:\t" + str_prec(a / 2.))
    print("2. / a:\t" + str_prec(2. / a))
    
    print("\n")
    
    print("-a:\t" + str_prec(-a))
    print("abs(c):\t" + str_prec(abs(c)))
    
    print("\n")
    
    print("pow(a,2.71):\t" + str_prec(pow(a,2.71)))
    print("a.sq()\t" + str_prec(a.sq()))
    print("a.sqrt()\t" + str_prec(a.sqrt()))
    print("a.cb()\t" + str_prec(a.cb()))
    print("a.cbrt()\t" + str_prec(a.cbrt()))
    print("a.exp()\t" + str_prec(a.exp()))
    print("a.log()\t" + str_prec(a.log()))
    
    print("a.sin()\t" + str_prec(a.sin()))
    print("a.cos()\t" + str_prec(a.cos()))
    print("a.tan()\t" + str_prec(a.tan()))
    # print("a.asin()\t" + str_prec(a.asin()))
    # print("a.acos()\t" + str_prec(a.acos()))
    # print("a.atan()\t" + str_prec(a.atan()))
    print("a.tanh()\t" + str_prec(a.tanh()))
    
    print("\n")
    print("\nTesting MCVectorData")
    print("\n------------------------\n")
    
    print("Manipulation\n")
    
    X = MCVectorData(np.array([2.3, 1.2, 0.7]), np.array([0.01, 0.01, 0.01]))
    Y = X+1.
    
    print("X:\n" + str_prec(X))
    print("Y:\n" + str_prec(Y))
    
    print("X + Y:\n" + str_prec(X+Y))
    print("X + 2.:\n" + str_prec(X+2.))
    print("2. + X:\n" + str_prec(2.+X))
    
    print("X + Y:\n" + str_prec(X+Y))
    print("X + 2.:\n" + str_prec(X+2.))
    print("2. + X:\n" + str_prec(2.+X))
    
    print("X / Y:\n" + str_prec(X/Y))
    print("X / 2.:\n" + str_prec(X/2.))
    print("2. / X:\n" + str_prec(2./X))
    
    print("X / Y:\n" + str_prec(X/Y))
    print("X / 2.:\n" + str_prec(X/2.))
    print("2. / X:\n" + str_prec(2./X))
    
    print("-X:\n" + str_prec(-X))
    print("abs(X):\n" + str_prec(X))
    
    print("pow(X,2.71):\n" + str_prec(pow(X,2.71)))
    print("X.sq():\n" + str_prec(X.sq()))
    print("X.sqrt():\n" + str_prec(X.sqrt()))
    print("X.cb():\n" + str_prec(X.cb()))
    print("X.cbrt():\n" + str_prec(X.cbrt()))
    print("X.exp():\n" + str_prec(X.exp()))
    print("X.log():\n" + str_prec(X.log()))
    
    print("X.sin():\n" + str_prec(X.sin()))
    print("X.cos():\n" + str_prec(X.cos()))
    print("X.tan():\n" + str_prec(X.tan()))
    # print("X.asin():\n" + str_prec(X.asin()))
    # print("X.acos():\n" + str_prec(X.acos()))
    # print("X.atan():\n" + str_prec(X.atan()))
    print("X.sinh():\n" + str_prec(X.sinh()))
    print("X.cosh():\n" + str_prec(X.cosh()))
    print("X.tanh():\n" + str_prec(X.tanh()))


if __name__ == '__main__':
    test_mcdata()