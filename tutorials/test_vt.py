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

import sys, os, subprocess, glob, platform


vistrailspath = '/root/vistrails/vistrails/run.py'
vtapp = [sys.executable, vistrailspath]
zipapp = 'unzip'
if platform.system()=='Windows':
  vtapp = [sys.executable, 'C:\\Program Files\\VisTrails\\vistrails\\run.py']
  zipapp = 'C:\\Program Files\\VisTrails\\unzip.exe',

if platform.system()=='Darwin':
  vtapp = ['/Applications/VisTrails/Vistrails.app/Contents/MacOS/vistrails']

# Find .vt files
vtfiles = glob.glob('[d-z]*/*.vt')
vtfiles.sort()

# Extract workflow tags from vt files
workflows = []
for vt in vtfiles:
    cmd = [zipapp, '-c', vt, 'vistrail']
    xmltrail = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]
    for line in xmltrail.splitlines():
        if line.find('key="__tag__"') == -1:
            continue
        tagstart = line.find('value="')+len('value="')
        tagend = line.find('"', tagstart+1)
        tag = line[tagstart:tagend]
        if tag != 'cannot prune':   # this seems to be some auto-generated tag
            workflows.append( (vt,tag) )
#            print os.path.basename(vt) + ':"' + tag + '"'

# Test all tagged workflows
logfile = open('vttest.log', 'w')
for workflow in workflows:
    (fn,tag) = workflow
    fn = os.path.join(os.getcwd(),fn)
    if not os.path.exists(fn):
        print(fn,'does not exist!')
    cmd = vtapp + ['-b', fn+':'+tag]
    print(cmd)
    print(fn + ':"' + tag + '" ', end=' ')
    logfile.write(str(cmd)+'\n')
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = process.communicate()
    logfile.write(fn+':'+tag+'\n')
    for k in output:
        logfile.write(k)
    logfile.write('===============================================================\n')
    print('returned',process.returncode)
