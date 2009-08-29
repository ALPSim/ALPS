# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Get ALPS at http://alps.comp-phys.org/
#
##############################################################################

from core.modules.vistrails_module import Module, ModuleError, NotCacheable
import core.bundles
import core.modules.basic_modules
import core.modules.module_registry

import alpscore
import plots


from plots import PlotFile

basic = core.modules.basic_modules

##############################################################################


class AlpsEvaluate(alpscore.SystemCommandLogged):
    def get_appname(self):
        return self.getInputFromPort('application')
    def compute(self):
        cmdlist = [self.get_appname()]
        for port_name in self.inputPorts:
           if port_name != 'file' and port_name != 'application' and self.hasInputFromPort(port_name):
             cmdlist += ['--'+str(port_name),str(self.getInputFromPort(port_name))]
        infile = self.getInputFromPort('file').name
        cmdlist += [infile]
        self.execute(cmdlist)
        outfiles = {}
        for port_name in self.outputPorts:
            outfiles[port_name] = PlotFile
            outfiles[port_name].name = infile.replace('.out.xml', '.plot.' + str(port_name) + '.xml')
            self.setResult(port_name,outfiles[port_name])
    _input_ports = [('file',[basic.File]),
                    ('application',[basic.File])]

class EvaluateFullDiagT(AlpsEvaluate):
    def get_appname(self):
      return alpscore._get_path('fulldiag_evaluate')
    _input_ports = [('T_MIN',[basic.Float]),
                    ('T_MAX',[basic.Float]),
                    ('DELTA_T',[basic.Float]),
                    ('application',[basic.File],True)]
    _output_ports = [('energy',[PlotFile]),
                    ('free_energy',[PlotFile]),
                    ('entropy',[PlotFile]),
                    ('specific_heat',[PlotFile]),
                    ('uniform_susceptibility',[PlotFile]),
                    ('magnetization',[PlotFile]),
                    ('application',[basic.File],True)]

class EvaluateFullDiagH(AlpsEvaluate):
    def get_appname(self):
      return alpscore._get_path('fulldiag_evaluate')
    _input_ports = [('H_MIN',[basic.Float]),
                    ('H_MAX',[basic.Float]),
                    ('DELTA_H',[basic.Float]),
                    ('application',[basic.File],True)]
    _output_ports = [('energy',[PlotFile]),
                    ('free_energy',[PlotFile]),
                    ('entropy',[PlotFile]),
                    ('specific_heat',[PlotFile]),
                    ('uniform_susceptibility',[PlotFile]),
                    ('magnetization',[PlotFile]),
                    ('application',[basic.File],True)]


def initialize(): pass

def selfRegister():

  reg = core.modules.module_registry.get_module_registry()

  reg.add_module(AlpsEvaluate,namespace="Evaluation",name="Evaluate",abstract=True)
  reg.add_module(EvaluateFullDiagT,namespace="Evaluation")
  reg.add_module(EvaluateFullDiagH,namespace="Evaluation")
  