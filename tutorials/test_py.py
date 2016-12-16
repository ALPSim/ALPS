from __future__ import print_function
# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Bela Bauer  <bauerb@gmail.com>
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

import sys, os, glob, subprocess

# The python interpreter running this script
pycmd = sys.executable

# Store current working directory
cwd = os.getcwd()

# Find .py files in subdirectories
pyfiles = glob.glob('[d-z]*/*.py')
pyfiles.sort()

# Test all Python tutorials
logfile = open('pytest.log', 'w')
for f in pyfiles:
    if f.find('build_lattice.py') >= 0: # This is no tutorial and needs cmd line arguments
        continue
    print(f + ':', end=' ')
    dir = os.path.dirname(f)
    fn = os.path.basename(f)
    os.chdir(os.path.join(cwd,dir))
    process = subprocess.Popen([pycmd,fn], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = process.communicate()
    logfile.write(f+':'+'\n')
    errors = []
    for k in output:
        logfile.write(k.decode(sys.stdout.encoding))
    logfile.write('===============================================================\n')
    print('returned',process.returncode)
