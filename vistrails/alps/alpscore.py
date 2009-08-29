# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Get ALPS at http://alps.comp-phys.org/
#
##############################################################################

from core.configuration import ConfigurationObject
from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.system import list2cmdline
import core.bundles
import core.modules.basic_modules
import core.modules.module_registry
import core.requirements
import os
import os.path
import tempfile

from packages.controlflow.list_module import ListOfElements

configuration = ConfigurationObject(path=(None, str))

basic = core.modules.basic_modules

binpath = ''

##############################################################################

def _get_path(binary_file):
    print 'get:', binpath, binary_file
    if binpath != '': 
        return os.path.join(binpath, binary_file)
    else:
        return binary_file

class SystemCommand(Module):
    def execute(self,cmdline):
        cmd = list2cmdline(cmdline)
        print cmd
        result = os.system(cmd)
        if result <> 0:
           raise ModuleError(self, 'Execution failed')

class SystemCommandLogged(Module):
    def execute(self,cmdline):
        logfile = self.interpreter.filePool.create_file(suffix='.log')
        cmdline += ['>&',logfile.name]
        print cmdline
        print "In execute"
        cmd = list2cmdline(cmdline)
        print cmd
        result = os.system(cmd)
        self.setResult('log_file', logfile)  
        if result <> 0:
           cmdline = ['open',logfile.name]
           cmd = list2cmdline(cmdline)
           os.system(cmd)
           raise ModuleError(self, 'Execution failed')
    _output_ports = [('log_file',[basic.File])]


class OpenHTML(NotCacheable, SystemCommand):
    """ open the file using the system open command """
    def compute(self):
        cmdlist = ['open', '-a', 'Safari']
        if self.hasInputFromPort('file'):
           cmdlist += [self.getInputFromPort('file').name]
        if self.hasInputFromPort('files'):
           cmdlist += self.getInputFromPort('files')
        self.execute(cmdlist)
    _input_ports = [('file', [basic.File]),
                    ('files', [ListOfElements])]


def initialize(): pass

def selfRegister():
  print "registering"
  reg = core.modules.module_registry.get_module_registry()
  
  reg.add_module(SystemCommand,namespace="Tools",abstract=True)
  reg.add_module(SystemCommandLogged,namespace="Tools",abstract=True)
  
  reg.add_module(OpenHTML,namespace="Tools")
