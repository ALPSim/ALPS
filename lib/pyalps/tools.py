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
import shutil
import tempfile
import subprocess
import platform
import sys
import glob
import numpy as np
import h5py

import pyalps.pytools # the C++ conversion functions
from load import loadBinningAnalysis, loadMeasurements,loadEigenstateMeasurements, loadSpectra, loadIterationMeasurements, loadObservableList, loadProperties
from hlist import deep_flatten, flatten, depth
from dict_intersect import dict_intersect
from dataset import DataSet
from plot_core import read_xml as readAlpsXMLPlot

def make_list(infiles):
    if type(infiles) == list:
      return infiles
    else:
      return [infiles]

def size(lst):
    try:
      return len(lst)
    except:
      return 1

def list2cmdline(lst):
    """ convert a list of arguments to a valid commandline """
    if platform.system() == 'Windows':
      return '"%s"' % subprocess.list2cmdline(lst)
    else:
      return subprocess.list2cmdline(lst)

def executeCommand(cmdline):
    """ execute the command given as list of arguments """
    cmd = list2cmdline(cmdline)
    print cmd
    return os.system(cmd)

def executeCommandLogged(cmdline,logfile):
    """ execute the command given as list of arguments and store the result into the log file """
    if platform.system() == 'Windows':
      cmdline += ['2>',logfile]
    else:
      cmdline += ['>&',logfile]
    return executeCommand(cmdline)

def runApplication(appname, parmfile, Tmin=None, Tmax=None, writexml=False, MPI=None, mpirun='mpirun'):
    """ run an ALPS application """
    cmdline = []
    if MPI != None:
        cmdline += [mpirun,'-np',str(MPI)]
    cmdline += [appname]
    if MPI != None:
        cmdline += ['--mpi','--Nmax','1']
    cmdline += [parmfile]
    if Tmin:
      cmdline += ['--Tmin',str(Tmin)]
    if Tmax:
      cmdline += ['--TMax',str(TMax)]
    if writexml:
      cmdline += ['--write-xml']
    return (executeCommand(cmdline),parmfile.replace('.in.xml','.out.xml'))
    
def runDMFT(infiles):
    """ run an DMFT application """
    appname='dmft'
    cmdline = [appname]
    cmdline += make_list(infiles)
    return (executeCommand(cmdline))
    
def evaluateLoop(infiles, appname='loop', write_xml=False):
    """ evaluate results of the looper QMC application """
    cmdline = [appname,'--evaluate']
    if write_xml:
      cmdline += ['--write_xml']
    cmdline += make_list(infiles)
    return executeCommand(cmdline)

def evaluateSpinMC(infiles, appname='spinmc_evaluate', write_xml=False):
    """ evaluate results of the Spin MC application """
    cmdline = [appname]
    if write_xml:
      cmdline += ['--write_xml']
    cmdline += make_list(infiles)
    return executeCommand(cmdline)

def evaluateQWL(infiles, appname='qwl_evaluate', DELTA_T=None, T_MIN=None, T_MAX=None):
    """ evaluate results of the quantum Wang-Landau application """
    cmdline = [appname]
    if DELTA_T:
      cmdline += ['--DELTA_T',str(DELTA_T)]
    if T_MIN:
      cmdline += ['--T_MIN',str(T_MIN)]
    if T_MAX:
      cmdline += ['--T_MAX',str(T_MAX)]
    cmdline += make_list(infiles)
    res = executeCommand(cmdline)
    if res != 0:
      raise Excpetion("Execution error in evaluateQWL: " + str(res))
    datasets = []
    for infile in infiles:
      datasets.append([])
      ofname = infile.replace('.out.xml', '.plot.*.xml')
      for fn in glob.glob(ofname):
        dataset = readAlpsXMLPlot(fn)
        datasets[-1].append(dataset)
        ylabel = dataset.props['ylabel']
    return datasets

def evaluateFulldiagVersusT(infiles, appname='fulldiag_evaluate', DELTA_T=None, T_MIN=None, T_MAX=None, H=None):
    """ evaluate results of the fulldiag application """
    cmdline = [appname]
    if DELTA_T != None:
      cmdline += ['--DELTA_T',str(DELTA_T)]
    if T_MIN != None:
      cmdline += ['--T_MIN',str(T_MIN)]
    if T_MAX != None:
      cmdline += ['--T_MAX',str(T_MAX)]
    if H != None:
      cmdline += ['--H',str(H)]
    cmdline += make_list(infiles)
    res = executeCommand(cmdline)
    if res != 0:
      raise Exception("Execution error in evaluateFulldiagVersusT: " + str(res))
    datasets = []
    for infile in infiles:
      datasets.append([])
      ofname = infile.replace('.out.xml', '.plot.*.xml')
      for fn in glob.glob(ofname):
        dataset = readAlpsXMLPlot(fn)
        datasets[-1].append(dataset)
        ylabel = dataset.props['ylabel']
    return datasets

def evaluateFulldiagVersusH(infiles, appname='fulldiag_evaluate', DELTA_H=None, H_MIN=None, H_MAX=None, T=None):
    """ evaluate results of the fulldiag application """
    cmdline = [appname,'--versus', 'h']
    if DELTA_H != None:
      cmdline += ['--DELTA_H',str(DELTA_H)]
    if H_MIN != None:
      cmdline += ['--H_MIN',str(H_MIN)]
    if H_MAX != None:
      cmdline += ['--H_MAX',str(H_MAX)]
    if T != None:
      cmdline += ['--T',str(T)]
    cmdline += make_list(infiles)
    res = executeCommand(cmdline)
    if res != 0:
      raise Exception("Execution error in evaluateFulldiagVersusH: " + str(res))
    datasets = []
    for infile in infiles:
      datasets.append([])
      ofname = infile.replace('.out.xml', '.plot.*.xml')
      for fn in glob.glob(ofname):
        dataset = readAlpsXMLPlot(fn)
        datasets[-1].append(dataset)
        ylabel = dataset.props['ylabel']
    return datasets

       
def inVistrails():
    """ returns True if called from within VisTrails """
    in_vistrails=True
    try:
      import core.modules.basic_modules
    except:
      in_vistrails=False
    return in_vistrails
    
def xslPath():
    """ return the path to the ALPS.xsl stylesheet """
    if inVistrails():
      if platform.system()=='Darwin':
        return os.path.join(sys.exec_prefix,'../Resources/lib/xml/ALPS.xsl')
      if platform.system()=='Windows':
        return os.path.join(sys.exec_prefix,'..','lib','xml','ALPS.xsl')
    return pyalps.pytools.search_xml_library_path("ALPS.xsl")
    
    
def copyStylesheet(dir):
    """ copy the ALPS.xsl stylesheet to the specified directory """
    target = os.path.join(dir,'ALPS.xsl')
    if not os.path.exists(target):
      shutil.copyfile(xslPath(), target)

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
    for fn in  glob.glob(filename.replace('.in.xml','.out.*')) + glob.glob(filename.replace('.in.xml','.clone*')):
      os.remove(fn)

def generateSeed():
    """ generate a random seed based on the current time
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
    if os.path.exists(fname+'.out.xml'):
      os.remove(fname+'.out.xml')
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
      taskname = base_name+'.task'+str(count)
      f.write('  <TASK status="new">\n')
      f.write('    <INPUT file="'+taskname+'.in.xml"/>\n')
      f.write('    <OUTPUT file="'+taskname+'.out.xml"/>\n')
      f.write('  </TASK>\n')
      writeTaskXMLFile(os.path.join(dirname,taskname+'.in.xml'),p)

    f.write('</JOB>\n')
    f.close()
    if (dirname==''):
      copyStylesheet('.')
    else:
      copyStylesheet(dirname)
    return fname+'.in.xml'
        

def writeParameterFile(fname,parms):
    """ This function writes a text input file for simple ALPS applications
    """
    f = file(fname,'w')
    for key in parms:
      value = parms[key]
      if type(value) == str:
        f.write(str(key)+' = "' + value + '"\n')
      else:
        f.write(str(key)+' = ' + str(value) + '\n')
    f.close()
    return fname


def recursiveGlob(dirname,pattern):
    ret = glob.glob(os.path.join(dirname, pattern))
    for d in os.listdir(dirname):
        if os.path.isdir(d):
            ret += recursiveGlob(os.path.join(dirname, d), pattern)
    return ret
    
def getResultFiles(dirname='.',pattern=None,prefix=None):
    """ get all result files matching the given pattern or prefix """
    if prefix!= None and pattern != None:
      raise Exception("Cannot define both prefix and pattern")
    if prefix == None: prefix = '*'
    if pattern == None:
      pattern = prefix+'.task*.out.xml'
      res=recursiveGlob(dirname, pattern)
      if len(res)==0:
        pattern = prefix+'.h5'
        res=recursiveGlob(dirname, pattern)
    else:
      res = recursiveGlob(dirname, pattern)
    return res

def groupSets(groups, for_each = []):
    dd = depth(groups)

    if dd > 1:
        hgroups = flatten(groups, -1)
        hgroups_idcs = hgroups.indices()
    else:
        hgroups = [groups]
        hgroups_idcs = [0]

    for idx in hgroups_idcs:
        sets = hgroups[idx]
    
        for_each_sets = {}
        for iset in sets:
            fe_par_set = tuple((iset.props[m] for m in for_each))
        
            if fe_par_set in for_each_sets:
                for_each_sets[fe_par_set].append(iset)
            else:
                for_each_sets[fe_par_set] = [iset]
    
        hgroups[idx] = for_each_sets.values()

    if dd > 1:
        return groups
    else:
        return hgroups[0]

def collectXY(sets,x,y,foreach=[]):
      foreach_sets = {}
      for iset in flatten(sets):
          if iset.props['observable'] != y:
              continue
          
          fe_par_set = tuple((iset.props[m] for m in foreach))
          
          if fe_par_set in foreach_sets:
              foreach_sets[fe_par_set].append(iset)
          else:
              foreach_sets[fe_par_set] = [iset]
      
      for k,v in foreach_sets.items():
          common_props = dict_intersect([q.props for q in v])
          res = DataSet()
          res.props = common_props
          for im in range(0,len(foreach)):
              m = foreach[im]
              res.props[m] = k[im]
          res.props['xlabel'] = x
          res.props['ylabel'] = y
          
          for data in v:
              if len(data.y)>1:
                  res.props['line'] = '.'
              xvalue = np.array([data.props[x] for i in range(len(data.y))])
              if len(res.x) > 0 and len(res.y) > 0:
                  res.x = np.concatenate((res.x, xvalue ))
                  res.y = np.concatenate((res.y, data.y))
              else:
                  res.x = xvalue
                  res.y = data.y
          
          order = np.argsort(res.x) #, kind = 'mergesort')
          res.x = res.x[order]
          res.y = res.y[order]
          res.props['label'] = ''
          for im in range(0,len(foreach)):
              res.props['label'] += '%s = %s ' % (foreach[im], k[im])
          
          foreach_sets[k] = res
      return foreach_sets.values()

def groupSets(groups, for_each = []):
    dd = depth(groups)

    if dd > 1:
        hgroups = flatten(groups, -1)
        hgroups_idcs = hgroups.indices()
    else:
        hgroups = [groups]
        hgroups_idcs = [0]

    for idx in hgroups_idcs:
        sets = hgroups[idx]

        for_each_sets = {}
        for iset in sets:
            fe_par_set = tuple((iset.props[m] for m in for_each))

            if fe_par_set in for_each_sets:
                for_each_sets[fe_par_set].append(iset)
            else:
                for_each_sets[fe_par_set] = [iset]

        hgroups[idx] = for_each_sets.values()

    if dd > 1:
        return groups
    else:
        return hgroups[0]

def subtract_spectrum(s1,s2,tolerance=1e-12):
    res = pyalps.DataSet()
    res.props = s1.props

    for i in range(len(s1.x)):
        remove = False
        for j in range(len(s2.x)):
            if abs(s1.x[i]-s2.x[j]) < tolerance and abs(s1.y[i]-s2.y[j]) < tolerance:
                remove = True
                break
        if not remove:
            res.x = np.append(res.x,s1.x[i])
            res.y = np.append(res.y,s1.y[i])
    
    return res

def save_parameters(filename, dict):
    f1 = h5py.File(filename, 'w')
    subgroup = f1.create_group('/parameters')
    
    for key in dict.keys():
        if(type(dict[key])==str):
            tid = h5py.h5t.C_S1.copy()
            tid.set_size(len(dict[key]))
        elif(type(dict[key])==int):
            tid = h5py.h5t.NATIVE_INT32.copy()
        else:
            tid = h5py.h5t.NATIVE_DOUBLE.copy()
        dset=subgroup.create_dataset(key, (), tid)
        dset[...] = dict[key]
    f1.close()
