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
import os
import os.path

from PyQt4 import QtCore, QtGui
from packages.spreadsheet.basic_widgets import SpreadsheetCell
from packages.spreadsheet.spreadsheet_cell import QCellWidget
import packages.spreadsheet


from packages.controlflow.list_module import ListOfElements

configuration = ConfigurationObject(path=(None, str))

basic = core.modules.basic_modules

binpath = ''

##############################################################################

def _get_path(binary_file):
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

class TextCell(SpreadsheetCell):
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
        self.display(TextCellWidget, (fileValue,))

class TextCellWidget(QCellWidget):
    """
  TextCellWidget has a QTextEdit to display HTML files
    
    """
    def __init__(self, parent=None):
        """ TextCellWidget(parent: QWidget) -> TextCellWidget
        Create a text cell without a toolbar and without editing capabilities
        
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

def selfRegister():
  print "registering"
  reg = core.modules.module_registry.get_module_registry()
  
  reg.add_module(SystemCommand,namespace="Tools",abstract=True)
  reg.add_module(SystemCommandLogged,namespace="Tools",abstract=True)
  
  reg.add_module(OpenHTML,namespace="Tools")

  reg.add_module(TextCell,namespace="Tools")
  reg.add_input_port(TextCell, "Location", packages.spreadsheet.basicWidgets.CellLocation)
  reg.add_input_port(TextCell, "File", basic.File)
