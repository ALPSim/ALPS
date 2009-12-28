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
import system
import glob

from core.modules.vistrails_module import ModuleError
from plots import PlotFile
from pyalps.plot_core import *
from dataset import DataSets, ResultFiles

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
        resultdir = basic.Directory
        resultdir.name = os.path.dirname(result.name)
        if self.hasInputFromPort('vtl'):
            filename = os.path.join(resultdir.name,'workflow.vtl')
            file_ = open(filename,'w')
            file_.write(self.getInputFromPort('vtl'))
            file_.close()
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
        self.setResult('output_file', result)
        self.setResult('output_dir', resultdir)
        
    _input_ports = [('input_file', [basic.File]),
                    ('tmin', [basic.Integer]),
                    ('tmax', [basic.Integer]),
                    ('continue', [basic.Boolean]),
                    ('application', [basic.File]),
                    ('vtl', [basic.String]),
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
           if port_name != 'files' and port_name != 'file' and port_name != 'application' and self.hasInputFromPort(port_name):
             cmdlist += ['--'+str(port_name),str(self.getInputFromPort(port_name))]
        rf = self.getInputFromPort('files')
        infiles = [x.props['filename'] for x in rf]
        print "Files", infiles
        cmdlist += infiles
        print cmdlist
        self.execute(cmdlist)
        datasetmap = {}
        datasets = []
        for infile in infiles:
          print "Looking at file", infile
          ofname = infile.replace('.out.xml', '.plot.*.xml')
          l = glob.glob(ofname)
          print "Have plots", l
          for fn in l:
            dataset = read_xml(fn)
            datasets.append(dataset)
            ylabel = dataset.props['ylabel']
            print ylabel
            if datasetmap.has_key(ylabel):
              datasetmap[ylabel].append(dataset)
            else:
              datasetmap[ylabel] = [dataset]
              
        for (port_name,ylabel) in self.plots:
          if datasetmap.has_key(ylabel):
            self.setResult(port_name,datasetmap[ylabel])
            print "For ",port_name,len(datasetmap[ylabel])
          else:
            self.setResult(port_name,[])
            print "Nothing for",port_name
          
    _input_ports = [('file',[basic.File],True),
                    ('files',[ResultFiles]),
                    ('application',[basic.File])]
    _output_ports = [('all',[DataSets])]
    appname = ''
    options = []

class EvaluateFullDiagT(AlpsEvaluate):
    appname = 'fulldiag_evaluate'
    _input_ports = [('T_MIN',[basic.Float]),
                    ('T_MAX',[basic.Float]),
                    ('DELTA_T',[basic.Float]),
                    ('application',[basic.File],True)]
    plots = [('energy','Energy Density'), 
             ('free_energy','Free Energy Density'),
             ('entropy','Entropy Density'),
             ('specific_heat','Specific Heat per Site'),
             ('uniform_susceptibility','Uniform Susceptibility per Site'),
             ('magnetization','Magnetization per Site'),
             ('compressibility','Compressibility per Site'),
             ('particle_number','Particle number per Site')]


class EvaluateFullDiagH(AlpsEvaluate):
    appname = 'fulldiag_evaluate'
    options = ['--versus', 'h']
    _input_ports = [('H_MIN',[basic.Float]),
                    ('H_MAX',[basic.Float]),
                    ('DELTA_H',[basic.Float]),
                    ('application',[basic.File],True)]
    plots = [('energy','Energy Density'), 
             ('free_energy','Free Energy Density'),
             ('entropy','Entropy Density'),
             ('specific_heat','Specific Heat per Site'),
             ('uniform_susceptibility','Uniform Susceptibility per Site'),
             ('magnetization','Magnetization per Site'),
             ('compressibility','Compressibility per Site'),
             ('particle_number','Particle number per Site')]


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
    plots = [('energy','Energy Density'), 
             ('free_energy','Free Energy Density'),
             ('entropy','Entropy Density'),
             ('specific_heat','Specific Heat per Site'),
             ('uniform_susceptibility','Uniform Susceptibility per Site'),
             ('staggered_structure_factor','Staggered Structure Factor per Site'),
             ('uniform_structure_factor','Uniform Structure Factor per Site')]



  
def initialize(): pass

def register_parameters(type, ns="Applications"):
  reg = core.modules.module_registry.get_module_registry()
  reg.add_module(type,namespace=ns)
  reg.add_output_port(type, "value", type)


def register_application(type):
  reg = core.modules.module_registry.get_module_registry()
  reg.add_module(type,namespace="Applications")
  reg.add_input_port(type,'application',[basic.File],False)

def register_evaluation(type):
  reg = core.modules.module_registry.get_module_registry()
  reg.add_module(type,namespace="Applications")
  reg.add_output_port(type,'all',[DataSets])
  for (port_name,ylabel) in type.plots:
    reg.add_output_port(type,port_name,[DataSets])
  reg.add_output_port(type,'log_file',[basic.File])
  
  
def selfRegister():

  reg = core.modules.module_registry.get_module_registry()
  
  register_parameters(system.SimulationID)
  
  reg.add_module(AlpsApplication,namespace="Applications")
  
  register_application(AppSpinMC)
  register_application(AppLoop)
  register_application(AppWorm)
  register_application(AppDirLoopSSE)
  register_application(AppFullDiag)
  register_application(AppSparseDiag)
  register_application(AppDMRG)
  register_application(AppQWL)
  
  reg.add_module(AlpsEvaluate,namespace="Applications",abstract=True)
  
  register_evaluation(EvaluateFullDiagT)
  register_evaluation(EvaluateFullDiagH)
  register_application(EvaluateLoop)
  register_evaluation(EvaluateQWL)
  
  reg.add_module(system.LatticeModel,namespace="Applications")
  reg.add_module(system.MonteCarloSimulation,namespace="Applications")
  reg.add_module(system.DiagonalizationSimulation,namespace="Applications")
  reg.add_module(system.DMRGSimulation,namespace="Applications")
