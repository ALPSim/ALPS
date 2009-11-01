# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Get ALPS at http://alps.comp-phys.org/
#
##############################################################################

from core.modules.vistrails_module import Module, ModuleError, NotCacheable
import core.bundles
import core.modules.basic_modules
import core.modules.module_registry
import os
import os.path
import tempfile
import copy
import glob

import parameters
import alpscore
import system
from parameters import Parameters
from packages.controlflow.list_module import ListOfElements

basic = core.modules.basic_modules

##############################################################################



class MakeParameterFile(Module):
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
             print type(res), res
             res.write(f);
         f.close()
         self.setResult('file', o)
         self.setResult('simulationid',os.path.basename(o.name))
     _input_ports = [('parms', [Parameters]),
                     ('simulationid',[system.SimulationID])]
     _output_ports=[('file', [basic.File]),
                    ('simulationid',[system.SimulationID])]


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
        l = glob.glob(name)
        self.setResult('value',l)
    def compute(self):
      self.expand(self.getInputFromPort('input_file').name)
    _input_ports = [('input_file',[basic.File])]
    _output_ports = [('value',[ListOfElements])]

class GetRunFiles(Module):
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
         self.setResult('value',glob.glob(os.path.join(dirname,prefix+ '.task' + tasks + '.out.run' +runs)))
     _input_ports = [('dir',[basic.Directory]), 
                     ('prefix',[basic.String]), 
                     ('tasks',[basic.String]),
                     ('runs',[basic.String])]
     _output_ports = [('value',[ListOfElements])]

class GetResultFiles(Module):
     def compute(self):
         prefix = '*'
         tasks = '*'
         d = self.getInputFromPort('dir')
         dirname = d.name
         if (self.hasInputFromPort('tasks')):
           tasks = self.getInputFromPort('tasks')
         self.setResult('value',glob.glob(os.path.join(dirname,prefix+ '.task' + tasks + '.out.xml')))
     _input_ports = [('dir',[basic.Directory]), 
                      ('prefix',[basic.String]), 
                      ('tasks',[basic.String])]
     _output_ports = [('value',[ListOfElements])]

class Convert2XML(alpscore.SystemCommandLogged):
    def compute(self):
        input_file = self.getInputFromPort('input_file')
        self.execute([alpscore._get_path('convert2xml')] + input_file)
        olist = []
        for q in input_file:
          olist.append(q + '.xml')
        self.setResult('value', olist)
    _input_ports = [('input_file', [ListOfElements])]
    _output_ports = [('value', [ListOfElements]),
                     ('log_file',[basic.File])]
 
class Convert2Text(alpscore.SystemCommand):
    def compute(self):
        input_file = self.getInputFromPort('input_file')
        output_file = self.interpreter.filePool.create_file(suffix='.txt')
        self.execute([alpscore._get_path('convert2xml'), input_file.name, '>' , output_file.name])
        self.setResult('output_file', output_file)
    _input_ports = [('input_file', [basic.File])]
    _output_ports = [('output_file', [basic.File])]

class XML2HTML(alpscore.SystemCommand):
    def compute(self):
        input_file = self.getInputFromPort('input_file')
        output_file = self.interpreter.filePool.create_file(suffix='.html')
        cmdlist = [alpscore._get_path('xslttransform')]
        if self.hasInputFromPort('stylesheet'):
          cmdlist += [self.getInputFromPort('stylesheet').name]
        cmdlist += [input_file.name, '>' , output_file.name]
        self.execute(cmdlist)
        self.setResult('output_file', output_file)
    _input_ports = [('input_file', [basic.File]),
                    ('stylesheet',[basic.File])]
    _output_ports = [('output_file', [basic.File])]


class PackSimulationResults(alpscore.SystemCommandLogged):
    def compute(self):
        o = self.interpreter.filePool.create_file(suffix='.tar.gz')
        if self.hasInputFromPort("dir"):
          dirname = self.getInputFromPort("dir").name
        self.execute(['cd', dirname,';', 'tar','czf', o.name, '*'])
        self.setResult("archive", o)
    _input_ports = [('dir', [basic.Directory])]
    _output_ports = [('archive', [basic.File]),
                     ('log_file',[basic.File])]


class GetSimName:
    def get_sim_name(self,dirname):
        l = glob.glob(os.path.join(dirname,'*.out.xml'))
        return l[0]


class UnpackSimulationResults(alpscore.SystemCommandLogged,GetSimName):
    def compute(self):
        o = self.interpreter.filePool.create_file()
        os.unlink(o.name)
        os.mkdir(o.name)
        dir = basic.Directory
        dir.name = o.name
        input_file = self.getInputFromPort("archive")
        self.execute(['cd', o.name,';', 'tar','xzf', input_file.name])
        o.name = self.get_sim_name(dir.name)
        self.setResult("output_file",o)
        self.setResult("output_dir",dir)
    _input_ports = [('archive', [basic.File])]
    _output_ports = [('output_file', [basic.File]),
                     ('output_dir', [basic.Directory]),
                     ('log_file',[basic.File])]        
        
        
class GetSimulationInDir(basic.Module,GetSimName):
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




def initialize(): pass

def selfRegister():

  reg = core.modules.module_registry.get_module_registry()
  
  reg.add_module(MakeParameterFile,namespace="Tools")
  reg.add_module(Parameter2XML,namespace="Tools")
  
  reg.add_module(Glob,namespace="Tools",abstract=True)
  reg.add_module(GetRunFiles,namespace="Tools")
  reg.add_module(GetResultFiles,namespace="Tools")
  
  reg.add_module(Convert2XML,namespace="Tools")
  reg.add_module(Convert2Text,namespace="Tools")
  reg.add_module(XML2HTML,namespace="Tools")

  reg.add_module(PackSimulationResults,namespace="Tools")
  reg.add_module(UnpackSimulationResults,namespace="Tools")

  reg.add_module(GetSimulationInDir,namespace="Tools")

  reg.add_module(PickFileFromList,namespace="Tools")
