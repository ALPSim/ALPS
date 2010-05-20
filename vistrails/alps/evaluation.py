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

from core.modules.vistrails_module import Module, ModuleError, NotCacheable
import core.bundles
import core.modules.basic_modules
import core.modules.module_registry

import alpscore
import plots
import tools
import os

from plots import PlotFile

basic = core.modules.basic_modules

##############################################################################


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


class EvaluateLoop(alpscore.SystemCommandLogged,tools.GetSimName):
    def compute(self):
        name = self.get_sim_name(self.getInputFromPort('dir').name)
        self.execute([alpscore._get_path('loop'),'--evaluate',name])
        self.setResult('dir',self.getInputFromPort('dir'))
    _input_ports = [('dir', [basic.Directory])]
    _output_ports = [('dir', [basic.Directory])]
    
class EvaluateSpinMC(alpscore.SystemCommandLogged,tools.GetSimName):
    def compute(self):
        name = self.get_sim_name(self.getInputFromPort('dir').name)
        self.execute([alpscore._get_path('spinmc'),name])
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

  reg.add_module(AlpsEvaluate,namespace="Evaluation",abstract=True)
  reg.add_module(EvaluateLoop,namespace="Evaluation",abstract=True)
  reg.add_module(EvaluateSpinMC,namespace="Evaluation",abstract=True)
  reg.add_module(EvaluateQWL,namespace="Evaluation",abstract=True)
  