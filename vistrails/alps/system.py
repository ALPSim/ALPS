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
import lattices 
import models

from PyQt4 import QtCore, QtGui
from packages.spreadsheet.basic_widgets import SpreadsheetCell
from packages.spreadsheet.spreadsheet_cell import QCellWidget
import packages.spreadsheet

basic = core.modules.basic_modules

##############################################################################

class SimulationID(basic.String):
   """ a simulation ID """

class LatticeModel(parameters.Parameters): 
   """ the simulation parameters, conistsing of model, lattice, and other parameters """
   def compute(self):
       self.parms = {}
       self.defaults = {}
       self.update_parms_from_port('lattice')
       self.update_parms_from_port('model')
       self.setOutput()
   _input_ports = [('lattice', [lattices.LatticeParameters]),
                    ('model', [models.ModelParameters])]
   _output_ports=[('value', [parameters.SystemParameters]),
                  ('value_as_string', [basic.String])]

class DiagonalizationSimulation(parameters.Parameters):
   """ a module collecting the typical input parameters for exact diagonalization """
   def compute(self):
       self.parms = {}
       self.defaults = {}
       for port_name in self.inputPorts:
          self.update_parms_from_port(port_name)
       self.setOutput()
   _input_ports = [('system', [parameters.SystemParameters]),
                    ('conserved', [parameters.ConservedQuantumnumbers]),
                    ('measurements',[parameters.CustomMeasurements])]
   _output_ports=[('value', [parameters.SystemParameters]),
                  ('value_as_string', [basic.String])]


class MonteCarloSimulation(parameters.Parameters):
   """ a module collecting the typical input parameters for a Monte Carlo simulation """
   def compute(self):
       self.parms = {}
       self.defaults = {}
       for port_name in self.inputPorts:
          self.update_parms_from_port(port_name)
       self.setOutput()
   _input_ports = [('system', [parameters.SystemParameters]),
                    ('mcparms', [parameters.MonteCarloParameters]),
                    ('temperature',[parameters.Temperature]),
                    ('measurements',[parameters.MonteCarloMeasurements])]
   _output_ports=[('value', [parameters.SystemParameters]),
                  ('value_as_string', [basic.String])]

class TextEditCell(SpreadsheetCell):
    """
    RichTextCell is a custom Module to view HTML files
    
    """
    def compute(self):
        """ compute() -> None
        Dispatch the HTML contents to the spreadsheet
        """
        if self.hasInputFromPort("File"):
            fileValue = self.getInputFromPort("File")
        else:
            fileValue = None
        self.display(TextEditCellWidget, (fileValue,))

class TextEditCellWidget(QCellWidget):
    """
    RichTextCellWidget has a QTextBrowser to display HTML files
    
    """
    def __init__(self, parent=None):
        """ RichTextCellWidget(parent: QWidget) -> RichTextCellWidget
        Create a rich text cell without a toolbar
        
        """
        QCellWidget.__init__(self, parent)
        self.setLayout(QtGui.QVBoxLayout(self))
        self.browser = QtGui.QTextEdit()
        self.layout().addWidget(self.browser)
        self.browser.setReadOnly(True)
 #       self.browser.controlBarType = None

    def updateContents(self, inputPorts):
        """ updateContents(inputPorts: tuple) -> None
        Updates the contents with a new changed in filename
        
        """
        (fileValue,) = inputPorts
        if fileValue:
            try:
                fi = open(fileValue.name, "r")
            except IOError:
                self.browser.setText("Cannot load the text file!")
                return            
            self.browser.setText(fi.read())
            fi.close()
        else:
            self.browser.setText("No text file is specified!")


def initialize(): pass

def register_parameters(type, ns="System"):
  reg = core.modules.module_registry.get_module_registry()
  reg.add_module(type,namespace=ns)
  reg.add_output_port(type, "value", type)
  reg.add_output_port(type, "value_as_string", basic.String)

def selfRegister():


  reg = core.modules.module_registry.get_module_registry()

  register_parameters(SimulationID)
  reg.add_module(LatticeModel,namespace="System")
  reg.add_module(MonteCarloSimulation,namespace="System")
  reg.add_module(DiagonalizationSimulation,namespace="System")

  reg.add_module(TextEditCell,namespace="System")
  reg.add_input_port(TextEditCell, "Location", packages.spreadsheet.basicWidgets.CellLocation)
  reg.add_input_port(TextEditCell, "File", basic.File)
