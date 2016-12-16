#! /usr/bin/env python
# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2009-2010 by Matthias Troyer <troyer@phys.ethz.ch> 
#                            Jan Gukelberger
#                            Adrian Feiguin
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

from __future__ import print_function
import sys

ns = sys.argv[1]
nedges = str(int(ns) - 1)

print('<LATTICES>')
print(' <GRAPH name = "open chain lattice with special edges" dimension="1" vertices="' + ns + '" edges="' + nedges + '">')
print('  <VERTEX id="1" type="0"><COORDINATE>0</COORDINATE></VERTEX>')
for i in range(2,int(ns)):
    print('  <VERTEX id="' + str(i) + '" type="1"><COORDINATE>' + str(i-1) + '</COORDINATE></VERTEX>')
print(' <VERTEX id="' + str(ns) + '" type="0"><COORDINATE>' + str(int(ns)-1) + '</COORDINATE></VERTEX>')
for i in range(1,int(ns)):
    nn = str(i + 1)
    print('  <EDGE source="' + str(i) + '" target="' + nn + '" id="' + str(i) + '" type="0" vector="1"/>')

print(' </GRAPH>')
print('</LATTICES>')
