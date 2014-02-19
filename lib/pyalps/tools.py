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
import math
import scipy.stats
import copy

import pyalps.hdf5 as h5

from pyalps.pytools import convert2xml, hdf5_name_encode, hdf5_name_decode, rng
import pyalps.pytools # the C++ conversion functions
from load import loadBinningAnalysis, loadMeasurements,loadEigenstateMeasurements, loadSpectra, loadIterationMeasurements, loadMPSIterations, loadObservableList, loadDMFTIterations, loadProperties, in_vistrails, log, Hdf5Loader
from hlist import deep_flatten, flatten, depth
from dict_intersect import dict_intersect
from dataset import DataSet
from floatwitherror import FloatWithError
from plot_core import read_xml as readAlpsXMLPlot
from dict_intersect import *
from natural_sort import natural_sort
import alea
import scipy.interpolate

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
      return subprocess.list2cmdline(lst)
    else:
      return subprocess.list2cmdline(lst)

def executeCommand(cmdline):
    """ execute the command given as list of arguments """
    cmd = list2cmdline(cmdline)
    log(cmd)
    # proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # sout, serr = proc.communicate() # serr should be empty
    # log(sout)
    # return proc.returncode
    return subprocess.call(cmd, shell=True)

def executeCommandLogged(cmdline,logfile):
    """ execute the command given as list of arguments and store the result into the log file """
    # subprocess is stupid: the interpretation of a list of args depends on shell=True|False
    # I'm using shell=True for backward compatibility to os.system
    cmd = list2cmdline(cmdline)
    log(cmd)
    return subprocess.call(cmd, shell=True, stdout=open(logfile, 'w'), stderr=subprocess.STDOUT)

def runApplication(appname, parmfiles, T=None, Tmin=None, Tmax=None, writexml=False, MPI=None, mpirun='mpirun'):
    """ run an ALPS application 
    
        This function runs an ALPS application. The parameers are:
        
        appname: the name of the application
        parmfile: the name of the main XML input file
        writexml: optional parameter, to be set to True if all results should be written to the XML files in addition to the HDF5 files
        T: time limit of MC simulation
        Tmin: optional parameter specifying the minimum time between checks whether a MC simulatio is finished
        Tmax: optional parameter specifying the maximum time between checks whether a MC simulatio is finished
        MPI: optional parameter specifying the number of processes to be used in an MPI simulation. MPI is not used if this parameter is left at ots default value  None.
        mpirun: optional parameter giving the name of the executable used to laucnh MPI applications. The default is 'mpirun'
    """
    if isinstance(parmfiles, str):
      parmfiles = [parmfiles];

    for parmfile in parmfiles:
      cmdline = []
      if MPI != None:
          cmdline += [mpirun,'-np',str(MPI)]
      cmdline += [appname]
      if MPI != None:
          cmdline += ['--mpi']
          if appname in ['sparsediag','fulldiag','dmrg']:
              cmdline += ['--Nmax','1']
      cmdline += [parmfile]
      if T:
        cmdline += ['-T',str(T)]
      if Tmin:
        cmdline += ['--Tmin',str(Tmin)]
      if Tmax:
        cmdline += ['--TMax',str(Tmax)]
      if writexml:
        cmdline += ['--write-xml']
      if parmfile.find('.xml') != -1:
        return (executeCommand(cmdline),parmfile.replace('.in.xml','.out.xml'))  # no iteration for xml i/o
      if parmfile.find('.h5') != -1:
        executeCommand(cmdline);
   

def runDMFT(infiles,apppath=''):
    """ run the ALPS DMFT application 
    
        The ALPS DMFT application does not (yet) use the standard ALPS input files and scheduler. Thus there is a separate function to call it. 
        This function takes one mandatory parameter: a single input file or a list of input files.
        Optional parameter apppath allows setting the path to the binary.
    """
    appname='dmft'
    return (executeCommand([apppath+appname] + make_list(infiles)))
    
def evaluateLoop(infiles, appname='loop', write_xml=False):
    """ evaluate results of the looper QMC application 
    
        this function calls the evaluate tool of the looper application. Additionally evaluated results are written back into the files. Besides a list of result files it takes one optional argument:
        
        write_xml: if this optional argument is set to True, the results will also bw written to the XML files
    """
    cmdline = [appname,'--evaluate']
    if write_xml:
      cmdline += ['--write-xml']
    cmdline += make_list(infiles)
    return executeCommand(cmdline)

def evaluateSpinMC(infiles, appname='spinmc_evaluate', write_xml=False):
    """ evaluate results of the spinmc application 
    
        this function calls the evaluate tool of the spinmc application. Additionally evaluated results are written back into the files. Besides a list of result files it takes one optional argument:
        
        write_xml: if this optional argument is set to True, the results will also bw written to the XML files
        
        
    """
    cmdline = [appname]
    if write_xml:
      cmdline += ['--write-xml']
    cmdline += make_list(infiles)
    return executeCommand(cmdline)

def evaluateQWL(infiles, appname='qwl_evaluate', DELTA_T=None, T_MIN=None, T_MAX=None):
    """ evaluate results of the quantum Wang-Landau application 
    
        this function calls the evaluate tool of the quantum Wang-Landau application. Besides a list of result files it takes the following arguments:
        T_MIN: the lower end of the temperature range for which quantities are evaluated
        T_MAX: the upper end of the temperature range for which quantities are evaluated
        DELTA_T: the temperature steps to be used between T_MIN and T_MAX
        
        This function returns a list of lists of DataSet objects, for the various properties evaluated for each of the input files.
    """
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
    """ evaluate results of the fulldiag application as a function of temperature
    
        this function calls the evaluate tool of the fulldiag application and evaluates several quantities as a function of temperature. Besides a list of result files it takes the following arguments:
        T_MIN: the lower end of the temperature range for which quantities are evaluated
        T_MAX: the upper end of the temperature range for which quantities are evaluated
        DELTA_T: the temperature steps to be used between T_MIN and T_MAX
        H: (optional) the magnetic field at which all data should be evaluated
        
        This function returns a list of lists of DataSet objects, for the various properties evaluated for each of the input files.
    """
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
    """ evaluate results of the fulldiag application as a function of magnetic field h
    
        this function calls the evaluate tool of the fulldiag application and evaluates several quantities as a function of magnetic field. Besides a list of result files it takes the following arguments:
        H_MIN: the lower end of the field range for which quantities are evaluated
        H_MAX: the upper end of the temperature range for which quantities are evaluated
        DELTA_H: the field steps to be used between H_MIN and H_MAX
        T: the temperature field at which all data should be evaluated
        
        This function returns a list of lists of DataSet objects, for the various properties evaluated for each of the input files.
    """
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
    #return in_vistrails
    return (sys.prefix.find('VisTrails') != -1)

    
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

def generateSeed():
    """ generate a random seed based on the current time
    """
    now = datetime.datetime.now()
    baseseed = now.microsecond+1000000*now.second+60000000*now.minute
    baseseed = ((baseseed << 10) | (baseseed >> 22));
    return baseseed

def writeInputH5Files(filename_,params_list):
  """ This function writes the H5 input files for ALPS (NGS)

      The parameters are:
      1. filename_   : the base file name of the H5 files that will be written
      2. params_list : a list of python dicts containing the simulation parameters

      Ping Nang MA 
  """
  input_files_ = [];
  for index in range(len(params_list)):
    if filename_.find('.in.h5') != -1:
      this_filename_ = filename_; 
    else:
      this_filename_ = filename_ + '.task' + str(index+1) + '.in.h5';
    input_files_.append(this_filename_);
    oar = pyalps.hdf5.archive(this_filename_,'w');
    for key in params_list[index].keys():
      oar['/parameters/' + key] = params_list[index][key]
    del oar;
  return input_files_;

def getInputH5Files(prefix='*',pattern='.task*.in.h5'):
    return glob.glob(prefix+pattern);


def getParameters(infiles_):
   """ This function extracts the parameters from the H5 input files for ALPS (NGS)

       1. infiles_ : either a list of filenames, or just one filename, containing NGS parameters object.
       
       Will be returned as a list of python dicts (parameters).

       Ping Nang MA
   """
   if isinstance(infiles_, str):
     infiles_ = [infiles_];

   params = [];
   for infile_ in infiles_:
     iar = pyalps.hdf5.archive(infile_);
     params_dict = {};
     for key in iar.list_children('/parameters'):
       params_dict[key] = iar['/parameters/' + key];
     params.append(params_dict);

   return params;


def writeInputFiles(fname,parms, baseseed=None):
    """ This function writes the XML input files for ALPS
    
         Parameters are:
         fname: the base file name of the XML files that will be written
         parms: a list of dicts containing the simulation parameters
         baseseed: optional parameter giving a random number seed from which seeds for the individual simulations will be calculated. The default value is taken from the current time.
    
         The function returns the name of the main XML input file
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
      if not 'SEED' in p:
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
    """ This function writes a text input file for simple ALPS applications like DMFT
    
        The arguments are:
        
          filename: the name of the parameter file to be written
          parms: the parameter dict
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

def input2output(infile):
    if   infile.find('.in.h5') != -1:
      outfile = infile.replace('.in.h5', '.out.h5');
    elif infile.find('.out.h5') != -1:
      outfile = infile;
    elif infile.find('.h5') != -1:
      outfile = infile.replace('.h5', '.out.h5');
    elif infile.find('.in.xml') != -1:
      outfile = infile.replace('.in.xml', '.out.xml');
    elif infile.find('.out.xml') != -1:
      outfile = infile;
    elif infile.find('.xml') != -1:
      outfile = infile.replace('.xml', '.out.xml');
    else:
      outfile = infile;
    return outfile;


def recursiveGlob(dirname,pattern):
    ret = glob.glob(os.path.join(dirname, pattern))
    for d in os.listdir(dirname):
        d = os.path.join(dirname, d)
        if os.path.isdir(d):
            ret += recursiveGlob(d, pattern)
    return ret

def getResultFiles(dirname='.',pattern=None,prefix=None,format=None):
    """ get all result files matching the given pattern or prefix 
    
        This function returns a list of all ALPS result files matching a given pattern, starting recursively from a given directory.
        The pattern can be either specificed by giving a prefix for the files, which is then augmented with the default ALPS file name suffixes. 
        ALternatively a fiull custom regular expression pattern can be specified.
        
        The paramters are:
        
        dirname: The directory from which to start the recursive search, defaulting to the current working directory.
        pattern: a regular expression pattern resricting the files to be matches
        prefix: a pattern which the start of the file names has to match. This will be augmented by the standard ALPS file name endings '*.task*.out.xml' or '*.h5' to form the full pattern.
        
        The function returns a list of filenames
    """
    if prefix!= None and pattern != None:
      raise Exception("Cannot define both prefix and pattern")
    if prefix == None: prefix = '*'
    if pattern == None:
      if format == None:
        pattern = prefix+'.task*.out.xml'
        res=recursiveGlob(dirname, pattern)
        if len(res)==0:
          pattern = prefix+'*.task*.out.h5'
          res=recursiveGlob(dirname, pattern)
          if len(res)==0:
            pattern = prefix+'*h5'
            res=recursiveGlob(dirname, pattern)
      else:
        if   format == 'xml':
          pattern = prefix+'.task*.out.xml'
          res=recursiveGlob(dirname, pattern)
        elif format == 'hdf5':
          pattern = prefix+'.task*.out.h5'
          res=recursiveGlob(dirname, pattern)
          if len(res)==0:
            pattern = prefix+'*h5'
            res=recursiveGlob(dirname, pattern)
    else:
      res = recursiveGlob(dirname, pattern)
    replicas=recursiveGlob(dirname, prefix+'*replica*h5')
    res += replicas
    return res

def loadTimeSeries(outfile, observable=None):
  if isinstance(outfile,str):
    if isinstance(observable,str):
      outfile = outfile.replace('.xml','.h5');

      ar = pyalps.hdf5.archive(outfile);
      return ar['/simulation/results'][observable]['timeseries']['data']

def getMeasurements(outfiles_, observable=None, includeLog=False):
  measurements = [];

  if isinstance(outfiles_,str):
    outfiles_ = [outfiles_];

  for outfile in outfiles_:
    outfile = outfile.replace('.xml','.h5');
    ar = pyalps.hdf5.archive(outfile);
    if isinstance(observable, str):
      measurements.append(ar['/simulation/results'][observable])
    
    del ar;

  return measurements;

def checkSteadyState(sets=None, outfile=None, observable=None, confidenceInterval=0.6827, includeLog=False):
  if sets != None:
    results = []
    for iset in flatten(sets):
      iset.props['checkSteadyState'] = checkSteadyState(outfile=iset.props['filename'], observable=iset.props['observable'], confidenceInterval=confidenceInterval, includeLog=True);
      results.append(iset); 
    return results 

  else:
    ts  = pyalps.loadTimeSeries(outfile, observable);  ### y
    N   = ts.size;
    idx = scipy.linspace(1, N, N);                     ### x

    beta1 = scipy.polyfit(idx, ts, 1)[0];              ### slope
    
    ts_std    = np.std(ts, ddof=1);                          ### unbiased estimate of standard deviation in y
    beta1_std = math.sqrt((12.*ts_std*ts_std)/(N * (N*N-1)));   ### unbiased estimate of standard deviation in slope

    z  = abs(beta1/beta1_std);
    z0 = scipy.stats.norm.ppf((1.-confidenceInterval) + 0.5*(confidenceInterval));   
    
    result = z < z0;

    if not includeLog:
      return {'value': result};
    else:
      return {'value': result, 'props':{ 'outfile': outfile, 'observable': observable}, 'statistics': {'beta1' : {'value' : beta1, 'std' : beta1_std}, 'confidenceInterval' : confidenceInterval, 'z' : z, 'z0' : z0}};

def checkConvergence(sets):
  results = []
  for iset in flatten(sets):
    ar = pyalps.hdf5.archive(iset.props['filename']) 
    if ar['/simulation/results'][iset.props['observable']]['mean']['error_convergence'] == 0:
      result = True;
    else:
      result = False;
    del ar;
    iset.props['checkConvergence'] = result
    results.append(iset);
  return results;

def sendmail(recipients, sender=None, message='', subject='', attachment=None):
  message = 'Automatic email message from ALPS.\nDo not reply.\n\n' + str(message);
  subject = 'Automatic email message from ALPS. ' + str(subject); 

  command = ['echo', message, '|', 'mail', '-s', subject];
  if sender != None:
    command += ['-r', sender];
  if attachment != None:
    command += ['-a', attachment]; 
  command += [recipients];

  return pyalps.executeCommand(command);

def extract(appname, source, target):
  cmdline = [appname];
  cmdline += [source];
  cmdline += [target];
  executeCommand(cmdline);
  return; 


def collectXY(sets,x,y,foreach=[],ignoreProperties=False):
      """ collects specified data from a list of DataSet objects
         
          this function is used to collect data from a list of DataSet objects, to prepare plots or evaluation. The parameters are:
    
            sets:    the list of datasets
            x:       the name of the property or measurement to be used as x-value of the collected results 
            y:       the name of the property or measurement to be used as y-value of the collected results 
            foreach: an optional list of properties used for grouping the results. A separate DataSet object is created for each unique set of values of the specified parameers.
            ignoreProperties: setting ignoreProperties=True prevents collectXY() from collecting properties.
            
          The function returns a list of DataSet objects.
      """
      foreach_sets = {}
      for iset in flatten(sets):
          if iset.props['observable'] != y and not y in iset.props:
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
              if data.props['observable'] == y:
                  if len(data.y)>1:
                      res.props['line'] = '.'
                  xvalue = np.array([data.props[x] for i in range(len(data.y))])
                  if len(res.x) > 0 and len(res.y) > 0:
                      res.x = np.concatenate((res.x, xvalue ))
                      res.y = np.concatenate((res.y, data.y))
                  else:
                      res.x = xvalue
                      res.y = data.y
              elif not ignoreProperties:
                  res.props['line'] = '.'
                  xvalue = np.array([ data.props[x] ])
                  if len(res.x) > 0 and len(res.y) > 0:
                      res.x = np.concatenate((res.x, xvalue ))
                      res.y = np.concatenate((res.y, np.array([ data.props[y] ])))
                  else:
                      res.x = xvalue
                      res.y = np.array([ data.props[y] ])
          
          order = np.argsort(res.x, kind = 'mergesort')
          res.x = res.x[order]
          res.y = res.y[order]
          res.props['label'] = ''
          for im in range(0,len(foreach)):
              res.props['label'] += '%s = %s ' % (foreach[im], k[im])
          
          foreach_sets[k] = res
      return foreach_sets.values()

def paramsAtFixedY(sets,x,y,fixedY,foreach=[]):
  XYs = collectXY(sets,x,y,foreach);

  params = [];
  for XY in XYs:
    param = XY.props;
    xlabel = param['xlabel'];
    del param['observable'];
    del param['hdf5_path'];
    del param['label'];
    del param['xlabel'];
    del param['ylabel'];

    x = XY.x;
    y = [y.mean for y in XY.y];

    f = scipy.interpolate.interp1d(y,x);

    if isinstance(fixedY, int) or isinstance(fixedY, float):
      fixedY = [fixedY];
    
    for this_fixedY in fixedY:
      xnew = f(this_fixedY);
      this_param = copy.deepcopy(param);
      this_param[xlabel] = float(xnew);
      params.append(this_param);

  return params;

def groupSets(groups, for_each = []):
    """ groups a list of DataSet objects into a list of lists
        
        this function groups a list of DataSet objects into a list of lists, according to the values of the properties given in the for_ech argument. DataSet objects with the same values of the properties given in for_each are grouped together.
        The parameters are:
          data: the data to be grouped
          for_each: the properties according to which the data is grouped
    """
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
    """ subtract one spectrum from another 
    
        this function takes two DataSet objects as input and returns a new DataSet that contains all (x,y) values that occur in the first but not in the second.
        The intended use is for plotting of spectra, where states that occur in several quantum number sectors should sometimes be shown only once.
    """
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

def save_parameters(filename, parms):
    """ saves parameters from a dict into an HDF5 file
    
        this function saves parameters from a Python dict into an ALPS HDF5 file.
        
        The arguments are:
        
          filename: the name of the HDF5 file
          parms: the parameter dict
    """
    f1=h5.archive(filename, 'w')
    for key in parms.keys():
        f1['/parameters/'+key] = parms[key]

def runTEBD(infileList):
    """ run a TEBD application """
    resList=[]
    appname='tebd'
    for infile in infileList:
        cmdline = [appname]
        cmdline += [infile]
        resList.append(executeCommand(cmdline))
    return resList

def stringListToList(inList):
    """ Convert a string which is of the form [...] with ... a collection of lists of floats and floats separated by commas
    into the list [...]."""
    outList=[]
    #sort out strings vs. floats
    if type(inList)==str:
        #remove whitespace and initial/final []
        dum=inList[1:len(inList)-1].replace(' ','')
        #find number of bracketed items (they come in pairs)
        numbrackets=dum.count('[')
        if numbrackets==0 :
            unbracketed=map(float,dum.replace('[','').replace(']','').replace(' ','').split(','))
            for q in unbracketed:
                outList.append([q])
        elif numbrackets>0:
            count=0
            for i in range(numbrackets):
                #find the first bracketed pair starting at count
                startInd=dum.find('[',count)
                finishInd=dum.find(']', count)
                if startInd>count:
                    unbracketed=map(float,(dum[count:startInd-1].replace('[','').replace(']','')\
                    .replace(' ','').split(',')))
                    for q in unbracketed:
                        outList.append([q])
                outList.append(map(float,dum[startInd:finishInd+1].replace('[','').\
                replace(']','').replace(' ','').split(',')))
                count=finishInd+2
            if len(dum)-count>0:
                unbracketed=map(float,dum[count:len(dum)].replace('[','').replace(']','')\
                .replace(' ','').split(','))
                for q in unbracketed:
                    outList.append([q])
    else:
        outList=inList
        for i in range(len(outList)):
            if type(outList[i])!=list:
                outList[i]=[outList[i]]
    return outList

def stringListToListstring(inList):
    """ Convert a string which is of the form [...] with ... a collection of lists of strings and strings separated by commas
    into the list [...]."""

    outList=[]
    #sort out strings vs. floats
    if type(inList)==str:
        #remove whitespace and initial/final []
        dum=inList[1:len(inList)-1].replace(' ','').replace("'",'')
        #find number of bracketed items (they come in pairs)
        numbrackets=dum.count('[')
        if numbrackets==0 :
            unbracketed=dum.replace('[','').replace(']','').replace(' ','').split(',')
            for q in unbracketed:
                outList.append([q.replace("'",'')])
        elif numbrackets>0:
            count=0
            for i in range(numbrackets):
                #find the first bracketed pair starting at count
                startInd=dum.find('[',count)
                finishInd=dum.find(']', count)
                if startInd>count:
                    unbracketed=(dum[count:startInd-1].replace('[','').replace(']','')\
                    .replace(' ','').split(','))
                    for q in unbracketed:
                        outList.append([q.replace("'",'')])
                outList.append(dum[startInd:finishInd+1].replace('[','').\
                replace(']','').replace(' ','').split(','))
                count=finishInd+2
            if len(dum)-count>0:
                unbracketed=dum[count:len(dum)].replace('[','').replace(']','')\
                .replace(' ','').split(',')
                for q in unbracketed:
                    outList.append([q.replace("'",'')])
    else:
        outList=inList
        for i in range(len(outList)):
            if type(outList[i])!=list:
                outList[i]=[outList[i]]
    return outList


def writeTEBDfiles(parmsList, fileName):
    counter=0
    nmlList=[]
    for parms in parmsList:
        counter+=1
        #Set up the systemSettings portion of the nameList file
        systemSettingsString='systemSize='+str(parms['L'])
        systemSettingsString=systemSettingsString+", Hamitype='"+str(parms['MODEL'])+"'"
        systemSettingsString=systemSettingsString+", initialState='"+str(parms['INITIAL_STATE'])+"'"
        if (not 'TAUS' in parms) :
            systemSettingsString+=', rtp=.false.'
        else:
            systemSettingsString+=', rtp=.true.'
        #check for conserved quantum numbers
        if 'CONSERVED_QUANTUMNUMBERS' in parms:
            if str(parms['MODEL'])=='spin':
                if parms['CONSERVED_QUANTUMNUMBERS']=='Sz_total':
                    systemSettingsString+=", qswitch=.true."
                    systemSettingsString+=", qType='Sz_total'"
                else:
                    raise Exception("Only Sz_total may be conserved for spin models!")
                    systemSettingsString+=", qswitch=.false."
                    systemSettingsString+=", qType='Sz_total'"
            else:
                if parms['CONSERVED_QUANTUMNUMBERS']=='N_total':
                    systemSettingsString+=", qswitch=.true."
                    systemSettingsString+=", qType='N_total'"
                else:
                    raise Exception("Only N_total may be conserved for particle models!")
                    systemSettingsString+=", qswitch=.false."
                    systemSettingsString+=", qType='N_total'"
        else:
            systemSettingsString+=", qswitch=.false."
            systemSettingsString+=", qType='Sz_total'"

        #check for value of conserved quantum number
        #Sz
        if 'Sz_total' in parms:
            if str(parms['MODEL'])=='spin':
                #Add L*S to this number to convert it to an integer appropriate for TEBD
                parms['Sz_total']=float(parms['L'])*float(parms['local_S'])-float(parms['Sz_total'])
                #convert to an integer-complain if it doesn't work
                if int(parms['Sz_total'])==float(parms['Sz_total']):
                    systemSettingsString+=", totQ="+str(int(parms['Sz_total']))
                else:
                    raise Exception("Invalid Sz_total encountered!")
                    systemSettingsString+=", totQ=0"
            else:
                raise Exception("Sz_total only conserved for spin models!")
                systemSettingsString+=", totQ=0"
        #N
        elif 'N_total' in parms:
            if str(parms['MODEL'])=='spin':
                raise Exception("N_total only conserved for particle models!")
                systemSettingsString+=", totQ=0"
            else:
                systemSettingsString+=", totQ="+str(int(parms['N_total']))
        else:
            if 'INITIAL_STATE' in parms:
                if parms['INITIAL_STATE']=='kink':
                    systemSettingsString+=", totQ=0"
                else:
                    raise Exception("Value of conserved quantity not specified!  Use N_total for particles or Sz_total for magnetization!")
            else :
                raise Exception("Value of conserved quantity not specified!  Use N_total for particles or Sz_total for magnetization!")

        #check for openmp threading
        if 'NUM_THREADS' in parms:
            systemSettingsString+=', numThr='+str(parms['NUM_THREADS'])
        else:
            systemSettingsString+=', numThr=1'
        #check for rtp chi cutoff
        if 'CHI_LIMIT' in parms:
            systemSettingsString+=', chiLimit='+str(parms['CHI_LIMIT'])
        else:
            systemSettingsString+=', chiLimit=100'
        #check for rtp truncation error cutoff
        if 'TRUNC_LIMIT' in parms:
            systemSettingsString+=', truncLimit=%30.15E' % (float(parms['TRUNC_LIMIT']))
        else:
            systemSettingsString+=', truncLimit=1.0E-12'
        #check for rtp truncation error cutoff
        if 'SIMID' in parms:
            systemSettingsString+=', simId='+str(parms['SIMID'])
        else:
            systemSettingsString+=', simId='+str(counter)
        #check for verbose switch
        if 'VERBOSE' in parms:
            systemSettingsString+=", print_switch=."+str(parms['VERBOSE'])+"."
        else:
            systemSettingsString+=",  print_switch=.false."
        systemSettingsString+='\n'

        #Output system settings
        nmlfileName=fileName+str(counter)+'.nml'
        nmlfile=open(nmlfileName,'w')
        nmlfile.write('&SystemSettings\n')
        nmlfile.write(systemSettingsString)
        nmlfile.write('&end\n\n')

        #find out which model
        mymodel=parms['MODEL']

        #spin and boson models have additional inputs
        if mymodel=='spin':
            nmlfile.write('&spinp\n')
            nmlfile.write('spin='+str(parms['local_S'])+'\n')
            nmlfile.write('&end\n\n')
        elif mymodel=='boson Hubbard':
            nmlfile.write('&bosonp\n')
            nmlfile.write('nmax='+str(parms['Nmax'])+'\n')
            nmlfile.write('&end\n\n')


        #if the ground state is the initial state, output itp parameters
        if parms['INITIAL_STATE']=='ground':
            #find out lengths of chi, trunc, and convCriteria, if they exist
            if 'ITP_CHIS' in parms:
                #convert from strings of the form '[float, float]' to
                # [float, float] etc. for use with vistrails modules
                if type(parms['ITP_CHIS'])==str:
                    itpChiList=map(int,(parms['ITP_CHIS'].replace('[','').replace(']','').replace(' ','').split(',')))
                else:
                    itpChiList=parms['ITP_CHIS']
            else:
                itpChiList=[50]
            if 'ITP_DTS' in parms:
                if type(parms['ITP_DTS'])==str:
                    itpDtList=map(float,(parms['ITP_DTS'].replace('[','').replace(']','').replace(' ','').split(',')))
                else:
                    itpDtList=parms['ITP_DTS']
            else:
                itpDtList=[0.01]
            if 'ITP_CONVS' in parms:
                if type(parms['ITP_CONVS'])==str:
                    itpConvList=map(float,(parms['ITP_CONVS'].replace('[','').replace(']','').replace(' ','').split(',')))
                else:
                    itpConvList=parms['ITP_CONVS']
            else:
                itpConvList=[1.0E-8]
            numItp=len(itpChiList)
            if len(itpChiList)==0:
                numItp=1
                itpChiList=[50]
            if len(itpDtList)==0:
                itpDtList=[0.01]
            if len(itpConvList)==0:
                itpConvList=[1.0E-8]
            if len(itpDtList)>numItp:
                itpChiList[numItp:len(itpDtList)-1]=itpChiList[numItp-1]
                numItp=len(itpDtList)
            elif len(itpDtList)<numItp:
                itpDtList[len(itpDtList):numItp-1]=itpDtList[len(itpDtList)-1]    
            if len(itpConvList)>numItp:
                itpChiList[numItp:len(itpDtList)-1]=itpChiList[numItp-1]
                itpDtList[numItp:len(itpConvList)-1]=itpDtList[numItp-1]
                numItp=len(itpConvList)
            elif len(itpConvList)<numItp:
                itpConvList[len(itpConvList):numItp-1]=itpConvList[len(itpConvList)-1]    
            itpfileName=fileName+str(counter)+'_itp.dat'
            #Write itp data to namelist file
            nmlfile.write('&ITPsettings\n')
            nmlfile.write('numITP='+str(numItp)+", itpfilename='"+itpfileName+"'\n")
            nmlfile.write('&end\n\n')

            #Write ITP data to itp file
            itpfile=open(itpfileName,'w')
            chiString=''
            for s in itpChiList:
                chiString+='%16i' %s
            chiString=chiString+'\n'
            itpfile.write(chiString)
            dtString=''
            for s in itpDtList:
                dtString+='%30.15E' %s
            dtString=dtString+'\n'
            itpfile.write(dtString)
            convString=''
            for s in itpConvList:
                convString+='%30.15E' %s
            convString=convString+'\n'
            itpfile.write(convString)
            itpfile.close()

            #Write ITP Hamiltonian parameters to namelist
            if mymodel=='spin':
                #Set up spin parameters
                myJz=0.0
                myJxy=0.0
                myH=0.0
                myGamma=0.0
                myD=0.0
                myK=0.0
                if 'J' in parms:
                    if parms['J']!='J0':
                        myJz=float(parms['J'])
                        myJxy=float(parms['J'])
                if 'Jz' in parms:
                    myJz=float(parms['Jz'])
                if 'Jxy' in parms:
                    myJxy=float(parms['Jxy'])
                if 'H' in parms:
                    myH=float(parms['H'])
                if 'Gamma' in parms:
                    myGamma=float(parms['Gamma'])
                if 'D' in parms:
                    myD=float(parms['D'])
                if 'K' in parms:
                    myK=float(parms['K'])
                nmlfile.write('&sp\n')
                itpnmlString='spinP%Jz='
                itpnmlString+='%30.15E'%(myJz)
                itpnmlString+=', spinP%Jxy='
                itpnmlString+='%30.15E'%(myJxy)
                itpnmlString+=', spinP%h='
                itpnmlString+='%30.15E'%(myH)
                itpnmlString+=', spinP%gam='
                itpnmlString+='%30.15E'%(myGamma)
                itpnmlString+=', spinP%d='
                itpnmlString+='%30.15E'%(myD)
                itpnmlString+=', spinP%k='
                itpnmlString+='%30.15E\n'%(myK)
                nmlfile.write(itpnmlString)
                nmlfile.write('&end\n\n')
            elif mymodel=='boson Hubbard':
                #Set up boson Hubbard parameters
                myT=1.0
                myU=0.0
                myV=0.0
                myMu=0.0
                if 't' in parms:
                    myT=float(parms['t'])
                if 'U' in parms:
                    myU=float(parms['U'])
                if 'V' in parms:
                    myV=float(parms['V'])
                if 'mu' in parms:
                    myMu=float(parms['mu'])
                nmlfile.write('&bp\n')
                itpnmlString='bosonP%mu='
                itpnmlString+='%30.15E'%(myMu)
                itpnmlString+=', bosonP%t='
                itpnmlString+='%30.15E'%(myT)
                itpnmlString+=', bosonP%V='
                itpnmlString+='%30.15E'%(myV)
                itpnmlString+=', bosonP%U='
                itpnmlString+='%30.15E\n'%(myU)
                nmlfile.write(itpnmlString)
                nmlfile.write('&end\n\n')
            elif mymodel=='hardcore boson':
                #Set up boson Hubbard parameters
                myT=1.0
                myV=0.0
                myMu=0.0
                if 't' in parms:
                    myT=float(parms['t'])
                if 'V' in parms:
                    myV=float(parms['V'])
                if 'mu' in parms:
                    myMu=float(parms['mu'])
                nmlfile.write('&hcbp\n')
                itpnmlString='hcbosonp%mu='
                itpnmlString+='%30.15E'%(myMu)
                itpnmlString+=', hcbosonp%t='
                itpnmlString+='%30.15E'%(myT)
                itpnmlString+=', hcbosonp%V='
                itpnmlString+='%30.15E\n'%(myV)
                nmlfile.write(itpnmlString)
                nmlfile.write('&end\n\n')
            elif mymodel=='fermion Hubbard':
                #Set up fermion Hubbard parameters
                myT=1.0
                myU=0.0
                myV=0.0
                myMu=0.0
                if 't' in parms:
                    myT=float(parms['t'])
                if 'U' in parms:
                    myU=float(parms['U'])
                if 'V' in parms:
                    myV=float(parms['V'])
                if 'mu' in parms:
                    myMu=float(parms['mu'])
                nmlfile.write('&fp\n')
                itpnmlString='fermiP%mu='
                itpnmlString+='%30.15E'%(myMu)
                itpnmlString+=', fermiP%t='
                itpnmlString+='%30.15E'%(myT)
                itpnmlString+=', fermiP%V='
                itpnmlString+='%30.15E'%(myV)
                itpnmlString+=', fermiP%U='
                itpnmlString+='%30.15E\n'%(myU)
                nmlfile.write(itpnmlString)
                nmlfile.write('&end\n\n')
            elif mymodel=='spinless fermions':
                #Set up spinless fermions parameters
                myT=1.0
                myV=0.0
                myMu=0.0
                if 't' in parms:
                    myT=float(parms['t'])
                if 'V' in parms:
                    myV=float(parms['V'])
                if 'mu' in parms:
                    myMu=float(parms['mu'])
                nmlfile.write('&sfp\n')
                itpnmlString='sfermiP%mu='
                itpnmlString+='%30.15E'%(myMu)
                itpnmlString+=', sfermiP%t='
                itpnmlString+='%30.15E'%(myT)
                itpnmlString+=', sfermiP%V='
                itpnmlString+='%30.15E\n'%(myV)
                nmlfile.write(itpnmlString)
                nmlfile.write('&end\n\n')

        #If rtp is desired, set up the RTP output files
        if (not 'TAUS' in parms) :
            numQuenches=0
        else:

            #get quench data
            if type(parms['TAUS'])==str:
                myTaus=map(float,(parms['TAUS'].replace('[','').replace(']','').replace(' ','').split(',')))
            else:
                myTaus=parms['TAUS']
            numQuenches=len(myTaus)

            if 'NUMSTEPS' in parms:
                if type(parms['NUMSTEPS'])==str:
                    myNumsteps=map(int,(parms['NUMSTEPS'].replace('[','').replace(']','').replace(' ','').split(',')))
                else:
                    myNumsteps=parms['NUMSTEPS']
            else:
                myNumsteps=[100]

            if 'STEPSFORSTORE' in parms:
                if type(parms['STEPSFORSTORE'])==str:
                    mySfs=map(int,(parms['STEPSFORSTORE'].replace('[','').replace(']','').replace(' ','').split(',')))
                else:
                    mySfs=parms['STEPSFORSTORE']
            else:
                mySfs=[1]


            numParams=[]
            
            #find the number of parameters quenched in each quench
            if 'POWS' in parms:
                myPows=stringListToList(parms['POWS'])
                #check if any of the pows are list-valued
                for q in myPows:
                    if type(q)==list:
                        numParams.append(len(q))
                    else:
                        numParams.append(1)
            else :
                myPows=[0.0]
                numQuenches=1


            if 'GIS' in parms:
                myGis=stringListToList(parms['GIS'])
            else :
                myGis=[0.0]
                numQuenches=1

            if 'GFS' in parms:
                myGfs=stringListToList(parms['GFS'])
            else :
                myGfs=[0.0]
                numQuenches=1

            if 'GS' in parms:
                myGs=stringListToListstring(parms['GS'])
            else :
                myGs=['t']
                numQuenches=1


            rtpfileName=fileName+str(counter)+'_rtp.dat'
            #Write rtp data to namelist file
            nmlfile.write('&RTPsettings\n')
            nmlfile.write('numQuenches='+str(numQuenches)+", rtpfilename='"+rtpfileName+"'\n")
            nmlfile.write('&end\n\n')
            #Write RTP data to rtp file
            rtpfile=open(rtpfileName,'w')
            tauString=''
            for s in myTaus:
                tauString+='%30.15E' %s
            tauString=tauString+'\n'
            rtpfile.write(tauString)
            nsString=''
            for s in myNumsteps:
                nsString+='%16i' %s
            nsString=nsString+'\n'
            rtpfile.write(nsString)
            sfsString=''
            for s in mySfs:
                sfsString+='%16i' %s
            sfsString=sfsString+'\n'
            rtpfile.write(sfsString)

            #write the nparams to file
            nParamsString=''
            for s in numParams:
                nParamsString+='%16i' %s
            nParamsString=nParamsString+'\n'
            rtpfile.write(nParamsString)
            powString=''
            for s in myPows:
                for q in s:
                    powString+='%30.15E' %q
                powString=powString+'\n'
            rtpfile.write(powString)
            giString=''
            for s in myGis:
                for q in s:
                    giString+='%30.15E' %q
                giString=giString+'\n'
            rtpfile.write(giString)
            gfString=''
            for s in myGfs:
                for q in s:
                    gfString+='%30.15E' %q
                gfString=gfString+'\n'
            rtpfile.write(gfString)
            gsString=''
            for s in myGs:
                for q in s:
                    gsString+=q.ljust(10)
                gsString=gsString+'\n'
            rtpfile.write(gsString)
            rtpfile.close()


            #Write RTP Hamiltonian parameters to namelist
            if mymodel=='spin':
                #Set up spin parameters
                myJz=0.0
                myJxy=0.0
                myH=0.0
                myGamma=0.0
                myD=0.0
                myK=0.0
                if 'J' in parms:
                    if parms['J']!='J0':
                        myJz=float(parms['J'])
                        myJxy=float(parms['J'])
                if 'Jz' in parms:
                    myJz=float(parms['Jz'])
                if 'Jxy' in parms:
                    myJxy=float(parms['Jxy'])
                if 'H' in parms:
                    myH=float(parms['H'])
                if 'Gamma' in parms:
                    myGamma=float(parms['Gamma'])
                if 'D' in parms:
                    myD=float(parms['D'])
                if 'K' in parms:
                    myK=float(parms['K'])
                nmlfile.write('&sp\n')
                itpnmlString='spinP%Jz='
                itpnmlString+='%30.15E'%(myJz)
                itpnmlString+=', spinP%Jxy='
                itpnmlString+='%30.15E'%(myJxy)
                itpnmlString+=', spinP%h='
                itpnmlString+='%30.15E'%(myH)
                itpnmlString+=', spinP%gam='
                itpnmlString+='%30.15E'%(myGamma)
                itpnmlString+=', spinP%d='
                itpnmlString+='%30.15E'%(myD)
                itpnmlString+=', spinP%k='
                itpnmlString+='%30.15E\n'%(myK)
                nmlfile.write(itpnmlString)
                nmlfile.write('&end\n\n')
            elif mymodel=='boson Hubbard':
                #Set up boson Hubbard parameters
                myT=1.0
                myU=0.0
                myV=0.0
                myMu=0.0
                if 't' in parms:
                    myT=float(parms['t'])
                if 'U' in parms:
                    myU=float(parms['U'])
                if 'V' in parms:
                    myV=float(parms['V'])
                if 'mu' in parms:
                    myMu=float(parms['mu'])
                nmlfile.write('&bp\n')
                itpnmlString='bosonP%mu='
                itpnmlString+='%30.15E'%(myMu)
                itpnmlString+=', bosonP%t='
                itpnmlString+='%30.15E'%(myT)
                itpnmlString+=', bosonP%V='
                itpnmlString+='%30.15E'%(myV)
                itpnmlString+=', bosonP%U='
                itpnmlString+='%30.15E\n'%(myU)
                nmlfile.write(itpnmlString)
                nmlfile.write('&end\n\n')
            elif mymodel=='hardcore boson':
                #Set up boson Hubbard parameters
                myT=1.0
                myV=0.0
                myMu=0.0
                if 't' in parms:
                    myT=float(parms['t'])
                if 'V' in parms:
                    myV=float(parms['V'])
                if 'mu' in parms:
                    myMu=float(parms['mu'])
                nmlfile.write('&hcbp\n')
                itpnmlString='hcbosonp%mu='
                itpnmlString+='%30.15E'%(myMu)
                itpnmlString+=', hcbosonp%t='
                itpnmlString+='%30.15E'%(myT)
                itpnmlString+=', hcbosonp%V='
                itpnmlString+='%30.15E\n'%(myV)
                nmlfile.write(itpnmlString)
                nmlfile.write('&end\n\n')
            elif mymodel=='fermion Hubbard':
                #Set up fermion Hubbard parameters
                myT=1.0
                myU=0.0
                myV=0.0
                myMu=0.0
                if 't' in parms:
                    myT=float(parms['t'])
                if 'U' in parms:
                    myU=float(parms['U'])
                if 'V' in parms:
                    myV=float(parms['V'])
                if 'mu' in parms:
                    myMu=float(parms['mu'])
                nmlfile.write('&fp\n')
                itpnmlString='fermiP%mu='
                itpnmlString+='%30.15E'%(myMu)
                itpnmlString+=', fermiP%t='
                itpnmlString+='%30.15E'%(myT)
                itpnmlString+=', fermiP%V='
                itpnmlString+='%30.15E'%(myV)
                itpnmlString+=', fermiP%U='
                itpnmlString+='%30.15E\n'%(myU)
                nmlfile.write(itpnmlString)
                nmlfile.write('&end\n\n')
            elif mymodel=='spinless fermions':
                #Set up spinless fermions parameters
                myT=1.0
                myV=0.0
                myMu=0.0
                if 't' in parms:
                    myT=float(parms['t'])
                if 'V' in parms:
                    myV=float(parms['V'])
                if 'mu' in parms:
                    myMu=float(parms['mu'])
                nmlfile.write('&sfp\n')
                itpnmlString='sfermiP%mu='
                itpnmlString+='%30.15E'%(myMu)
                itpnmlString+=', sfermiP%t='
                itpnmlString+='%30.15E'%(myT)
                itpnmlString+=', sfermiP%V='
                itpnmlString+='%30.15E\n'%(myV)
                nmlfile.write(itpnmlString)
                nmlfile.write('&end\n\n')
        nmlList.append(nmlfileName)
    return nmlList



def select(inp,condition):
    data_ = []
    for ds in pyalps.flatten(inp):
        if condition(ds):
            data_.append(ds)
    return data_

def select_by_property(data,proplist):
    for k,v in proplist.items():
        data = select(data, lambda ds: ds.props[k]==v)
    return data

def values(data, key):    
    vals = []
    if type(key) == list:
        for ds in pyalps.flatten(data):
            keyv = tuple([ds.props[kk] for kk in key])
            if keyv not in vals:
                vals.append(keyv)
        return vals
    else:
        for ds in pyalps.flatten(data):
            if ds.props[key] not in vals:
                vals.append(ds.props[key])
        return np.sort(vals)

def mergeDataSets(dsets):
    props = dict_intersect([d.props for d in dsets])
    merged = copy.deepcopy(dsets.pop())
    for d in dsets:
        if (d.x != merged.x).any():   raise ValueError('cannot merge datasets: x values mismatch')
        if len(d.y) != len(merged.y): raise ValueError('cannot merge datasets: y shape mismatch')
        if isinstance(merged.y,alea.MCVectorData):
            merged.y.merge(d.y)
        else:
            for i in range(len(merged.y)):
                merged.y[i].merge(d.y[i])
    merged.props = props
    return merged

def mergeMeasurements(measurements):
    byname = {}
    for mset in measurements:
        for m in mset:
            key = m.props['observable']
            if key not in byname:   byname[key] = [m]
            else:                   byname[key].append(m)
    merged = [mergeDataSets(v) for v in byname.values()]
    return merged

def mergeMeasurementsFromFiles(files,respath='/simulation/realizations/0/clones/0/measurements'):
    ll = Hdf5Loader()
    meas = ll.ReadMeasurementFromFile(files,respath=respath)
    return mergeMeasurements(meas)

def saveMeasurements(measurements,outfile,respath='/simulation/results'):
    for m in measurements:
        path = respath+'/'+m.props['observable']
        if isinstance(m.y,alea.MCVectorData):
            m.y.save(outfile,path)
        elif isinstance(m.y,np.ndarray) and isinstance(m.y[0],alea.MCScalarData):
            m.y[0].save(outfile,path)
        elif isinstance(m.y,FloatWithError):
            h5f = h5.archive(fn, 'w')
            h5f[path+'/mean/value'] = np.array(m.y.mean)
            h5f[path+'/mean/error'] = np.array(m.y.error)
            try:
                h5f[path+'/jackknife'] = np.array(m.jacks)
            except AttributeError:
                pass

def SetLabels (data, proplist):
    """
    Set labels according to the properties given in 'proplist'.
    """
    if type(proplist) == str:
        proplist = [proplist]
    for d in flatten(data):
        if 'label' in d.props:
            del d.props['label']
        label_items = []
        for propname in proplist:
            try:
                label_items.append(propname+' = %s' % (d.props[propname]))
            except:
                pass
        d.props['label'] = ', '.join(label_items)
    return data

def CycleColors (data, foreach,
                 colors=['k','b','g','m','c','y']):
    """
    Cyclically assign colors to the lines/markers that will be used to
    display the DataSets, based on the properties in 'foreach'. This means
    that DataSet instances that have the same values for the properties
    in 'foreach' will receive the same color.
    """
    all = {}
    for q in flatten(data):
        key = tuple([q.props[k] for k in foreach])
        all[key] = ''

    icolor = 0
    for k in all.keys():
        all[k] = colors[icolor]
        icolor = (icolor+1)%len(colors)

    for q in flatten(data):
        key = tuple([q.props[k] for k in foreach])
        q.props['color'] = all[key]

    return data

def CycleMarkers (data, foreach,
                  markers=['s', 'o', '^', '>', 'v', '<', 'd', 'p', 'h', '+', 'x']):
    """
    Cyclically assign markers to the lines/markers that will be used to
    display the DataSets, based on the properties in 'foreach'. This means
    that DataSet instances that have the same values for the properties
    in 'foreach' will receive the same marker.
    """
    all = {}
    for q in flatten(data):
        key = tuple([q.props[k] for k in foreach])
        all[key] = ''

    imarker = 0
    for k in all.keys():
        all[k] = markers[imarker]
        imarker = (imarker+1)%len(markers)

    for q in flatten(data):
        key = tuple([q.props[k] for k in foreach])
        q.props['marker'] = all[key]
        if 'line' in q.props:
            ll = list(q.props['line'])
            ll[0] = all[key]
            q.props['line'] = ''.join(ll)
        else:
            q.props['line'] = all[key] + '-'

    return data

