# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Get ALPS at http://alps.comp-phys.org/
#
##############################################################################

import core.bundles
import core.modules.basic_modules
import core.modules.module_registry

import alpscore
import os

from core.modules.vistrails_module import ModuleError

basic = core.modules.basic_modules

##############################################################################


class AlpsApplication(alpscore.SystemCommandLogged):
    """ Runs an ALPS application for a given parameter file """

    def get_path(self,appname):
        if self.hasInputFromPort('num_processes') and  self.getInputFromPort('num_processes') > 1:
            return alpscore._get_path(appname+'_mpi')
        else: 
            return alpscore._get_path(appname)

    def get_app_name(self):
        if self.hasInputFromPort('application'):
          fn = self.getInputFromPort('application')
          return self.get_path(fn.name)
        else:
          if self.appname != '':
            return self.get_path(self.appname)
          else: 
             raise ModuleError(self, 'No application specified')

             
    def getoptions(self):
        options = []
        if self.hasInputFromPort('tmin'):
            options += ['--Tmin', str(self.getInputFromPort('tmin'))]
        if self.hasInputFromPort('tmax'):
            options += ['--Tmax', str(self.getInputFromPort('tmax'))]
        return options
    
    def compute(self):
        input_file = self.getInputFromPort('input_file')
        result = basic.File()
        result.name = input_file.name.replace('.in.xml', '.out.xml')
        an = self.get_app_name()
        if not os.path.isfile(an):
            raise ModuleError(self, "Application '%s' not existent" % an)
        cmdline = [an] + self.getoptions()
        if self.hasInputFromPort('continue'):
            cmdline += [result.name]
        else:
            cmdline += [input_file.name]
        print cmdline
        self.execute(cmdline)
        self.setResult('output_file', result)
        
    _input_ports = [('input_file', [basic.File]),
                    ('tmin', [basic.Integer]),
                    ('tmax', [basic.Integer]),
                    ('continue', [basic.Boolean]),
                    ('application', [basic.File]),
                    ('num_processes',[basic.Integer])
                    ]
    _output_ports = [('output_file', [basic.File]),
                     ('log_file',[basic.File])]
    appname=''
                         
class AppSpinMC(AlpsApplication):
    """Runs spinmc for given parameter file """
    appname = 'spinmc'

class SpinMC(AlpsApplication):
    """Runs spinmc for given parameter file """
    appname = 'spinmc'

class AppLoop(AlpsApplication):
    """Runs loop for given parameter file """
    def getoptions(self):
        options = []
        if self.hasInputFromPort('tmin'):
            options += ['--report-interval', str(self.getInputFromPort('tmin'))]
        return options
    appname = 'loop'

class AppDirLoopSSE(AlpsApplication):
    """Runs dirloop_sse for given parameter file """
    appname = 'dirloop_sse'

class AppWorm(AlpsApplication):
    """Runs worm for given parameter file """
    appname = 'worm'

class AppFullDiag(AlpsApplication):
    """Runs fulldiag for given parameter file """
    appname = 'fulldiag'

class AppSparseDiag(AlpsApplication):
    """Runs sparsediag for given parameter file """
    appname = 'sparsediag'

class AppDMRG(AlpsApplication):
    """Runs dmrg for given parameter file """
    appname = 'dmrg'

class AppQWL(AlpsApplication):
    """Runs qwl for given parameter file """
    appname = 'qwl'

  
def initialize(): pass

def selfRegister():

  reg = core.modules.module_registry.get_module_registry()
  
  reg.add_module(AlpsApplication,namespace="Applications")
  reg.add_module(AppSpinMC,namespace="Applications")
  reg.add_module(AppLoop,namespace="Applications")
  reg.add_module(AppWorm,namespace="Applications")
  reg.add_module(AppDirLoopSSE,namespace="Applications")
  reg.add_module(AppFullDiag,namespace="Applications")
  reg.add_module(AppSparseDiag,namespace="Applications")
  reg.add_module(AppDMRG,namespace="Applications")
  reg.add_module(AppQWL,namespace="Applications")
