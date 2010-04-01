# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Bela Bauer <bauerb@phys.ethz.ch>
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

import os.path
import datetime

def writeTaskXMLFile(filename,parms):
    f = file(filename,'w')
    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write('<?xml-stylesheet type="text/xsl" href="ALPS.xsl"?>\n')
    f.write('<SIMULATION xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://xml.comp-phys.org/2003/8/QMCXML.xsd">\n')
    f.write('  <PARAMETERS>\n')
    for key in parms:
      f.write('<PARAMETER name="'+str(key)+'">'+str(parms[key])+'</PARAMETER>\n')
    f.write('  </PARAMETERS>\n')
    f.write('</SIMULATION>\n')
    f.close()

def generateSeed():
    """ This function generates a random seed based on the current time
    """
    now = datetime.datetime.now()
    baseseed = now.microsecond+1000000*now.second+60000000*now.minute
    baseseed = ((baseseed << 10) | (baseseed >> 22));
    return baseseed

def writeInputFiles(fname,parms, baseseed=None):
    """ This function writes the XML input files for ALPS
    """
    dirname = os.path.dirname(fname)
    base_name = os.path.basename(fname)
    f = file(fname+'.in.xml','w')
    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write('<?xml-stylesheet type="text/xsl" href="ALPS.xsl"?>\n')
    f.write('<JOB xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://xml.comp-phys.org/2003/8/job.xsd">\n')
    f.write('  <OUTPUT file="'+base_name+'.out.xml"/>\n')
   
    bits = 31;
    n = len(parms)
    while n>0:
      n /= 2
      bits -= 1

    if baseseed == None:
      baseseed = generateSeed()
     
    count = 0
    for p in parms:
      count += 1
      if not p.has_key('SEED'):
        seed = baseseed
        for j in range(0,32/bits+1):
          seed ^= ((count-1) << (j * bits))
          seed &= ((1<<30) | ((1<<30)-1))
          p['SEED'] = seed
      taskname = fname+'.task'+str(count)
      f.write('  <TASK status="new">\n')
      f.write('    <INPUT file="'+taskname+'.in.xml"/>\n')
      f.write('    <OUTPUT file="'+taskname+'.out.xml"/>\n')
      f.write('  </TASK>\n')
      writeTaskXMLFile(taskname+'.in.xml',p)

    f.write('</JOB>\n')
    f.close()
    
    return fname+'.in.xml'
    #alpscore.copy_stylesheet(dir.name)