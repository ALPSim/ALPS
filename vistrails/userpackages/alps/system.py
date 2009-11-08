# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Get ALPS at http://alps.comp-phys.org/
#
##############################################################################

from core.modules.vistrails_module import Module, ModuleError, NotCacheable
import core.bundles
import core.modules.basic_modules
import core.modules.module_registry

import parameters
import alpsparameters
import lattices 
import models

from PyQt4 import QtCore, QtGui
from packages.spreadsheet.basic_widgets import SpreadsheetCell
from packages.spreadsheet.spreadsheet_cell import QCellWidget
import packages.spreadsheet

from alpsparameters import SystemParameters
basic = core.modules.basic_modules

##############################################################################

class SimulationID(basic.String):
    """ a simulation ID """

class LatticeModel(parameters.Parameters): 
    """ the simulation parameters, conistsing of model, lattice, and other parameters """
    def compute(self):
        res=self.updateFromPort('parms',parameters.ParametersData({}))
        res=self.updateFromPort('lattice',res)
        res=self.updateFromPort('model',res)
        self.setOutput(res)
    _input_ports = [('lattice', [lattices.LatticeParameters]),
                     ('model', [models.ModelParameters])]
    _output_ports=[('value', [SystemParameters])]

class DiagonalizationSimulation(parameters.Parameters):
    """ a module collecting the typical input parameters for exact diagonalization """
    def compute(self):
        res = parameters.ParametersData({})
        for port_name in self.inputPorts:
           res=self.updateFromPort(port_name,res)
        self.setOutput(res)
    _input_ports = [('system', [SystemParameters]),
                     ('conserved', [alpsparameters.ConservedQuantumnumbers]),
                     ('measurements',[alpsparameters.CustomMeasurements])]
    _output_ports=[('value', [SystemParameters])]


class MonteCarloSimulation(parameters.Parameters):
    """ a module collecting the typical input parameters for a Monte Carlo simulation """
    def compute(self):
        res = parameters.ParametersData({})
        for port_name in self.inputPorts:
           res=self.updateFromPort(port_name,res)
        self.setOutput(res)
    _input_ports = [('system', [SystemParameters]),
                     ('mcparms', [alpsparameters.MonteCarloParameters]),
                     ('temperature',[alpsparameters.Temperature]),
                     ('measurements',[alpsparameters.MonteCarloMeasurements])]
    _output_ports=[('value', [SystemParameters])]


def initialize(): pass

def register_parameters(type, ns="System"):
  reg = core.modules.module_registry.get_module_registry()
  reg.add_module(type,namespace=ns)
  reg.add_output_port(type, "value", type)

def selfRegister():


  reg = core.modules.module_registry.get_module_registry()

  register_parameters(SimulationID)
  reg.add_module(LatticeModel,namespace="System")
  reg.add_module(MonteCarloSimulation,namespace="System")
  reg.add_module(DiagonalizationSimulation,namespace="System")

