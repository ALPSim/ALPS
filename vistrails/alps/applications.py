# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Copyright (C) 2009 - 2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
#                              Synge Todo <wistaria@comp-phys.org>
#
# Distributed under the Boost Software License, Version 1.0. (See accompany-
# ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
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

from packages.vtlcreator.init import VtlFileCreator
from core.modules.vistrails_module import ModuleError
from plots import PlotFile
from pyalps.plot_core import *
from dataset import DataSets, ResultFiles

basic = core.modules.basic_modules

##############################################################################


class RunAlpsApplication(alpscore.SystemCommandLogged):
    """ Runs an ALPS application for a given parameter file """

    def get_path(self,appname):
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

    def num_procs(self):
        if self.hasInputFromPort('num_processes'):
          np = self.getInputFromPort('num_processes')
        else:
          np = alpscore._get_default_mpi_procs()
        return np

    def getoptions(self):
        options = self.options
        if self.hasInputFromPort('tmin'):
            options += ['--Tmin', str(self.getInputFromPort('tmin'))]
        if self.hasInputFromPort('tmax'):
            options += ['--Tmax', str(self.getInputFromPort('tmax'))]
        if self.hasInputFromPort('write_xml'):
          if self.getInputFromPort('write_xml'):
            options += ['--write-xml']
        return options
    
    def compute(self):
        input_file = self.getInputFromPort('input_file')
        result = basic.File()
        result.name = input_file.name.replace('.in.xml', '.out.xml')
        resultdir = basic.Directory()
        resultdir.name = os.path.dirname(result.name)
        an = self.get_app_name()
        if not os.path.isfile(an):
            raise ModuleError(self, "Application '%s' not existent" % an)
        if self.num_procs()>1:
            cmdline = alpscore._get_mpi_run()+[str(self.num_procs()),an,"--mpi"]
        else:
            cmdline = [an]
        cmdline += self.getoptions()
        if self.hasInputFromPort('continue'):
            cmdline += [result.name]
        else:
            cmdline += [input_file.name]
        f = file(os.path.join(resultdir.name,'workflow.vtl'),'w')
        f.write(VtlFileCreator.generate_vtl(self.moduleInfo['locator'],self.moduleInfo['version'],self.moduleInfo['pipeline']))
        f.close()
        self.execute(cmdline)
        self.setResult('output_file', result)
        self.setResult('output_dir', resultdir)
        
    _input_ports = [('input_file', [basic.File]),
                    ('tmin', [basic.Integer]),
                    ('tmax', [basic.Integer]),
                    ('continue', [basic.Boolean]),
                    ('application', [basic.File]),
                    ('num_processes',[basic.Integer]),
                    ('write_xml',[basic.Boolean])
                    ]
    _output_ports = [('output_file', [basic.File]),
                     ('output_dir', [basic.Directory]),
                     ('log_file',[basic.File])]
    appname=''
    options=[]
                         
class RunSpinMC(RunAlpsApplication):
    """Runs spinmc for given parameter file """
    appname = 'spinmc'

class RunLoop(RunAlpsApplication):
    """Runs loop for given parameter file """
    def getoptions(self):
        options = ['--auto-evaluate']
        if self.hasInputFromPort('tmin'):
            options += ['--report-interval', str(self.getInputFromPort('tmin'))]
        return options
    appname = 'loop'

class RunDirLoopSSE(RunAlpsApplication):
    """Runs dirloop_sse for given parameter file """
    appname = 'dirloop_sse'

class RunWorm(RunAlpsApplication):
    """Runs worm for given parameter file """
    appname = 'worm'

class RunFullDiag(RunAlpsApplication):
    """Runs fulldiag for given parameter file """
    appname = 'fulldiag'
    options = ['--Nmax', '1']

class RunSparseDiag(RunAlpsApplication):
    """Runs sparsediag for given parameter file """
    appname = 'sparsediag'
    options = ['--Nmax', '1']

class RunDMRG(RunAlpsApplication):
    """Runs dmrg for given parameter file """
    appname = 'dmrg'
    options = ['--Nmax', '1']

class RunQWL(RunAlpsApplication):
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

class EvaluateFullDiagVersusT(AlpsEvaluate):
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


class EvaluateFullDiagVersusH(AlpsEvaluate):
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
  
#  register_parameters(system.SimulationID)
  
  reg.add_module(RunAlpsApplication,namespace="Applications")
  
  register_application(RunSpinMC)
  register_application(RunLoop)
  register_application(RunWorm)
  register_application(RunDirLoopSSE)
  register_application(RunFullDiag)
  register_application(RunSparseDiag)
  register_application(RunDMRG)
  register_application(RunQWL)
  
  reg.add_module(AlpsEvaluate,namespace="Applications",abstract=True)
  
  register_evaluation(EvaluateFullDiagVersusT)
  register_evaluation(EvaluateFullDiagVersusH)
  register_application(EvaluateLoop)
  register_evaluation(EvaluateQWL)
  
  register_parameters(system.SimulationName)

  reg.add_module(system.LatticeModel,namespace="Applications")
  reg.add_module(system.PrepareMonteCarlo,namespace="Applications")
  reg.add_module(system.PrepareDiagonalization,namespace="Applications")
  reg.add_module(system.PrepareDMRG,namespace="Applications")
