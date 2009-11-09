# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Get ALPS at http://alps.comp-phys.org/
#
##############################################################################

import core.bundles
import core.modules.basic_modules
import core.modules.module_registry

import alpscore
import plots
import tools
import os

from core.modules.vistrails_module import ModuleError
from plots import PlotFile

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
        if self.hasInputFromPort('num_processes') and  self.getInputFromPort('num_processes') > 1:
            cmdline = alspcore.config.mpirun+[self.getInputFromPort('num_processes')]+cmdline
            print "Using MPI", cmdline
        if self.hasInputFromPort('continue'):
            cmdline += [result.name]
        else:
            cmdline += [input_file.name]
        print cmdline
        self.execute(cmdline)
        resultdir = basic.Directory
        resultdir.name = os.path.dirname(result.name)
        self.setResult('output_file', result)
        self.setResult('output_dir', resultdir)
        
    _input_ports = [('input_file', [basic.File]),
                    ('tmin', [basic.Integer]),
                    ('tmax', [basic.Integer]),
                    ('continue', [basic.Boolean]),
                    ('application', [basic.File]),
                    ('num_processes',[basic.Integer])
                    ]
    _output_ports = [('output_file', [basic.File]),
    				 ('output_dir', [basic.Directory]),
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


class AlpsEvaluate(alpscore.SystemCommandLogged):
    def get_appname(self):
        if self.hasInputFromPort('application'):
          return alpscore._get_path(self.getInputFromPort('application'))
        else:
          if self.appname != '':
            return alpscore._get_path(self.appname)
          else: 
             raise ModuleError(self, 'No application specified')
          
    def compute(self):
        an = self.get_appname()
        if not os.path.isfile(an):
            raise ModuleError(self, "Application '%s' not existent" % an)
        cmdlist = [an]
        cmdlist += self.options
        for port_name in self.inputPorts:
           if port_name != 'file' and port_name != 'application' and self.hasInputFromPort(port_name):
             cmdlist += ['--'+str(port_name),str(self.getInputFromPort(port_name))]
        infile = self.getInputFromPort('file').name
        cmdlist += [infile]
        self.execute(cmdlist)
        outfiles = {}
        for port_name in self.outputPorts:
            of = PlotFile()
            of.name = infile.replace('.out.xml', '.plot.' + str(port_name) + '.xml')
            self.setResult(port_name,of)
    _input_ports = [('file',[basic.File]),
                    ('application',[basic.File])]
    appname = ''
    options = []

class EvaluateFullDiagT(AlpsEvaluate):
    appname = 'fulldiag_evaluate'
    _input_ports = [('T_MIN',[basic.Float]),
                    ('T_MAX',[basic.Float]),
                    ('DELTA_T',[basic.Float]),
                    ('application',[basic.File],True)]
    _output_ports = [('energy',[PlotFile]),
                    ('free_energy',[PlotFile]),
                    ('entropy',[PlotFile]),
                    ('specific_heat',[PlotFile]),
                    ('uniform_susceptibility',[PlotFile]),
                    ('magnetization',[PlotFile])]


class EvaluateFullDiagH(AlpsEvaluate):
    appname = 'fulldiag_evaluate'
    options = ['--versus', 'h']
    _input_ports = [('H_MIN',[basic.Float]),
                    ('H_MAX',[basic.Float]),
                    ('DELTA_H',[basic.Float]),
                    ('application',[basic.File],True)]
    _output_ports = [('energy',[PlotFile]),
                    ('free_energy',[PlotFile]),
                    ('entropy',[PlotFile]),
                    ('specific_heat',[PlotFile]),
                    ('uniform_susceptibility',[PlotFile]),
                    ('magnetization',[PlotFile])]


class EvaluateLoop(alpscore.SystemCommandLogged,tools.GetSimName):
    def compute(self):
        name = self.get_sim_name(self.getInputFromPort('dir').name)
        self.execute([alpscore._get_path('loop'),'--evaluate',name])
        self.setResult('dir',self.getInputFromPort('dir'))
    _input_ports = [('dir', [basic.Directory])]
    _output_ports = [('dir', [basic.Directory])]

class EvaluateQWL(AlpsEvaluate):
    appname = 'qwl_evaluate'
    _input_ports = [('T_MIN',[basic.Float]),
                    ('T_MAX',[basic.Float]),
                    ('DELTA_T',[basic.Float]),
                    ('application',[basic.File],True)]
    _output_ports = [('energy',[PlotFile]),
                    ('free_energy',[PlotFile]),
                    ('entropy',[PlotFile]),
                    ('specific_heat',[PlotFile]),
                    ('uniform_susceptibility',[PlotFile]),
                    ('staggered_structure_factor',[PlotFile]),
                    ('uniform_structure_factor',[PlotFile])]



  
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
  reg.add_module(AlpsEvaluate,namespace="Applications",abstract=True)
  reg.add_module(EvaluateFullDiagT,namespace="Applications")
  reg.add_module(EvaluateFullDiagH,namespace="Applications")
  reg.add_module(EvaluateLoop,namespace="Applications")
  reg.add_module(EvaluateQWL,namespace="Applications")
