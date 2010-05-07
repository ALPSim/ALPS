# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Copyright (C) 2009 - 2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
#                              Synge Todo <wistaria@comp-phys.org>,
#                              Bela Bauer <bauerb@phys.ethz.ch>
#
# Distributed under the Boost Software License, Version 1.0. (See accompany-
# ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
#
##############################################################################

from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.configuration import ConfigurationObject
import core.bundles
import core.modules.basic_modules
import core.modules.module_registry
import os
import os.path
import tempfile
import copy
import glob
import zipfile
import datetime

import platform
import codecs

import parameters
import alpscore
import system
from parameters import Parameters
from packages.controlflow.list_module import ListOfElements

import pyalps
import pyalps.pytools # the C++ conversion functions


basic = core.modules.basic_modules

##############################################################################

from pyalps import ResultFile
from dataset import ResultFiles
from dataset.dataset_exceptions import EmptyInputPort, InvalidInput

class UnzipDirectory(Module):
    _input_ports = [('zipfile',[basic.File])]
    _output_ports = [('output_dir', [basic.Directory])]
    
    def compute(self):
        o = self.interpreter.filePool.create_file()
        os.unlink(o.name)
        os.mkdir(o.name)
        dir = basic.Directory
        dir.name = o.name
        os.chdir(dir.name)
        
        input_file = self.getInputFromPort('zipfile').name
        zf = zipfile.ZipFile(input_file,'r')
        # ugly, but necessary in Python 2.5
        filelist = zf.namelist()
        for f in filelist:
            if f[-1] == '/':
                os.mkdir(f)
            else:
                open(f,'w').write(zf.read(f))
        # This will work when Vistrails moves to Python 2.6
        # zf.extractall()
        
        self.setResult('output_dir',dir)

class WriteParameterFile(Module):
     """Creates a parameter file.
     """
     def compute(self):
         o = self.interpreter.filePool.create_file()
         if self.hasInputFromPort('simulationid'):
             o = self.interpreter.filePool.create_file()
             o.name = os.path.join(os.path.dirname(o.name),self.getInputFromPort('simulationid'))
         f = file(o.name,'w')
         if self.hasInputFromPort('parms'):
           input_values = self.forceGetInputListFromPort('parms')
           for p in input_values:
             res = parameters.make_parameter_data(p)
             res.write(f);
         f.close()
         self.setResult('file', o)
         self.setResult('simulationid',os.path.basename(o.name))
     _input_ports = [('parms', [Parameters]),
                     ('simulationid',[basic.String])]
     _output_ports=[('file', [basic.File]),
                    ('simulationid',[basic.String])]


class WriteInputFiles(Module):
     """ This module writes the XML input files for ALPS
     """
         
     def compute(self):
         of = self.interpreter.filePool.create_file()
         os.unlink(of.name)
         os.mkdir(of.name)
         dir = basic.Directory
         dir.name = of.name
         o = self.interpreter.filePool.create_file()
         if self.hasInputFromPort('simulationid'):
           base_name = self.getInputFromPort('simulationid')
         else:
           base_name = os.path.basename(o.name)

         ofile = basic.File()
         ofile.name = os.path.join(dir.name,base_name + '.in.xml')

         if self.hasInputFromPort('parms'):
           input_values = self.forceGetInputListFromPort('parms')
           l = []
           for p in input_values:
             if isinstance(p,list):
               l += p
             else:
               l += [p]

           if self.hasInputFromPort('baseseed'):
             baseseed = self.getInputFromPort('baseseed')
           else:
             baseseed = pyalps.generateSeed()
             
           Module.annotate(self,{'baseseed':baseseed})
           
           pyalps.writeInputFiles(os.path.join(dir.name,base_name),l,baseseed)
           pyalps.copyStylesheet(dir.name)
           
         self.setResult("output_dir", dir)
         self.setResult("output_file", ofile)
     _input_ports = [('parms', [Parameters]),
                     ('baseseed',[basic.Integer],True),
                     ('simulationid',[basic.String])]
     _output_ports = [('output_file', [basic.File]),
                     ('output_dir', [basic.Directory])]

class Parameter2XML(alpscore.SystemCommandLogged):
    def compute(self):
        o = self.interpreter.filePool.create_file()
        os.unlink(o.name)
        os.mkdir(o.name)
        input_file = self.getInputFromPort("file")
        base_name = os.path.basename(input_file.name)
        dir = basic.Directory
        dir.name = o.name
        self.execute(['cd',o.name,';', alpscore._get_path('parameter2xml'),
                            input_file.name, base_name])
        # Things would be easier on our side if ALPS somehow
        # advertised all the files it creates. Right now, it will be
        # hard to make sure all temporary files that are created will
        # be cleaned up. We are trying to keep all the temporaries
        # inside a temp directory, but the lack of control over where
        # the files are created hurt.
        ofile = basic.File()
        ofile.name = os.path.join(o.name,base_name + '.in.xml')
        self.setResult("output_dir", dir)
        self.setResult("output_file", ofile)
    _input_ports = [('file', [basic.File]),
                    ('output_dir', [basic.File],True)]
    _output_ports = [('output_file', [basic.File]),
                     ('output_dir', [basic.Directory]),
                     ('log_file',[basic.File])]


class Glob(Module):
    def expand(self,name):
        l = recursive_glob(name)
        self.setResult('value',l)
    def compute(self):
      self.expand(self.getInputFromPort('input_file').name)
    _input_ports = [('input_file',[basic.File])]
    _output_ports = [('value',[ListOfElements])]

class GetCloneFiles(Module):
     def compute(self):
         tasks = '*'
         runs = '*[0-9]'
         prefix = '*'
         d = self.getInputFromPort('dir')
         dirname = d.name
         if (self.hasInputFromPort('tasks')):
           tasks = str(self.getInputFromPort('tasks'))
         if (self.hasInputFromPort('runs')):
           tasks = str(self.getInputFromPort('runs'))
         if (self.hasInputFromPort('prefix')):
           tasks = str(self.getInputFromPort('prefix'))+'*'
         files = glob.glob(os.path.join(dirname,prefix+ '.task' + tasks + '.out.run' +runs))
# remove HDF5 files from the list
         res = []
         for f in files:
           if f[-3:]!='.h5':
             res.append(f)
         self.setResult('value',res)
     _input_ports = [('dir',[basic.Directory]), 
                     ('prefix',[basic.String]), 
                     ('tasks',[basic.String]),
                     ('runs',[basic.String])]
     _output_ports = [('value',[ListOfElements])]

class GetResultFiles(Module):
    def compute(self):
        if self.hasInputFromPort('pattern') and self.hasInputFromPort('dir'):
            dir = self.getInputFromPort('dir').name
            # try several ways to match the pattern
            pattern = self.getInputFromPort('pattern')
            pattern = os.path.expanduser(pattern)
            pattern = os.path.expandvars(pattern)
            # result = glob.glob(pattern)
            result = pyalps.recursiveGlob(dir,pattern)
            self.setResult('value',result)
            self.setResult('resultfiles', [ResultFile(x) for x in result])
        else:
            prefix = '*'
            tasks = '*'
            d = self.getInputFromPort('dir')
            dirname = d.name
            if self.hasInputFromPort('tasks'):
                tasks = self.getInputFromPort('tasks')
            if self.hasInputFromPort('prefix'):
                prefix = self.getInputFromPort('prefix')
            result = pyalps.recursiveGlob(dirname, prefix+ '.task' + tasks + '.out.xml')
            self.setResult('value', result)
            self.setResult('resultfiles', [ResultFile(x) for x in result])

    _input_ports = [('dir',[basic.Directory]), 
        ('prefix',[basic.String]), 
        ('tasks',[basic.String]),
        ('pattern',[basic.String])]
    _output_ports = [('value',[ListOfElements]),
        ('resultfiles',[ResultFiles])]

class Convert2XML(Module):
    def compute(self):
        input_files = self.getInputFromPort('input_file')
        olist = []
        for f in input_files:
          print "Converting", f
          olist.append(pyalps.pytools.convert2xml(str(f)))
          pyalps.copyStylesheet(os.path.dirname(f))
        self.setResult('value', olist)
    _input_ports = [('input_file', [ListOfElements])]
    _output_ports = [('value', [ListOfElements])]
 
class Convert2Text(alpscore.SystemCommand):
    def compute(self):
        input_file = self.getInputFromPort('input_file')
        output_file = self.interpreter.filePool.create_file(suffix='.txt')
        self.execute([alpscore._get_path('convert2text'), input_file.name, '>' , output_file.name])
        self.setResult('output_file', output_file)
    _input_ports = [('input_file', [basic.File])]
    _output_ports = [('output_file', [basic.File])]

class ConvertXML2HTML(alpscore.SystemCommand):
    def compute(self):
        input_file = self.getInputFromPort('input_file')
        output_file = self.interpreter.filePool.create_file(suffix='.html')
        if platform.system() == 'Windows':
          cmdlist = ['msxsl.exe',input_file.name]
          if self.hasInputFromPort('stylesheet'):
            cmdlist += [self.getInputFromPort('stylesheet').name]
          else:
            cmdlist += [pyalps.xslPath()]
          cmdlist += ['-o', output_file.name]
        if platform.system() != 'Windows':
          cmdlist = [alpscore._get_path('xslttransform')]
          if self.hasInputFromPort('stylesheet'):
            cmdlist += [self.getInputFromPort('stylesheet').name]
          cmdlist += [input_file.name, '>' , output_file.name]
        self.execute(cmdlist)
        if platform.system() == 'Windows': # need to convert to UTF-8
          fin = codecs.open(output_file.name,"r","utf-16")
          u = fin.read()
          fin.close()
          fout = file(output_file.name,"w")
          fout.write(u.encode("utf-8"))
          fout.close()
        self.setResult('output_file', output_file)
    _input_ports = [('input_file', [basic.File]),
                    ('stylesheet',[basic.File])]
    _output_ports = [('output_file', [basic.File])]


class GetSimName:
    def get_sim_name(self,dirname):
        l = glob.glob(os.path.join(dirname,'*.out.xml'))
        return l[0]


class GetJobFile(basic.Module,GetSimName):
    def compute(self):
        dir = self.getInputFromPort("dir")
        o = basic.File
        o.name = self.get_sim_name(dir.name)
        self.setResult("file",o)
        self.setResult("output_file",o)
    _input_ports = [('dir', [basic.Directory])]
    _output_ports = [('file', [basic.File]),
                     ('output_file', [basic.File],True)]
                     
class PickFileFromList(basic.Module):
    def compute(self):
        f=basic.File()
        ind = 0
        if self.hasInputFromPort('index'):
          ind = self.getInputFromPort('index')
        f.name = self.getInputFromPort('files')[ind]
        self.setResult('file',f)
    _input_ports = [('files', [ListOfElements]),
                    ('index', [basic.Integer])]
    _output_ports = [('file', [basic.File])]

class ArchiveDirectory(basic.Module):
    _input_ports = [
        ('archive', [basic.String]),
        ('path', [basic.Directory])
    ]
    _output_ports = [('output', [basic.Directory])]
    
    def compute(self):
        replace_dict = {}
        if alpscore.config.check('archives'):
            replace_dict = basic.Dictionary.translate_to_python(alpscore.config.archives)
        else:
            raise InvalidInput("Check configuration!")
        
        if self.hasInputFromPort('archive'):
            path = self.getInputFromPort('archive')
            for k,v in replace_dict.items():
                for qq in v:
                    testpath = path.replace(k, qq)
                    print 'Testing',testpath
                    if os.path.exists(testpath):
                        dir = basic.Directory()
                        dir.name = testpath
                        self.setResult('output', dir)
                        return
        
        if self.hasInputFromPort('path'):
            lp = self.getInputFromPort('path').name
            if os.path.exists(lp):
                dir = basic.Directory()
                dir.name = lp
                self.setResult('output', dir)
                return
        
        raise InvalidInput("Can't locate file at any location.")

def initialize(): pass

def selfRegister():

  reg = core.modules.module_registry.get_module_registry()
  
  reg.add_module(WriteParameterFile,namespace="Tools",abstract=True)
  reg.add_module(Parameter2XML,namespace="Tools",abstract=True)
  reg.add_module(WriteInputFiles,namespace="Tools")
  
  reg.add_module(Glob,namespace="Tools",abstract=True)
  reg.add_module(GetCloneFiles,namespace="Tools")
  reg.add_module(GetResultFiles,namespace="Tools")
  
  reg.add_module(Convert2XML,namespace="Tools")
  reg.add_module(Convert2Text,namespace="Tools")
  reg.add_module(ConvertXML2HTML,namespace="Tools")

  reg.add_module(UnzipDirectory,namespace="Tools")
  reg.add_module(ArchiveDirectory,namespace="Tools")

  reg.add_module(GetJobFile,namespace="Tools")

  reg.add_module(PickFileFromList,namespace="Tools")
