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

from core.configuration import ConfigurationObject
from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.system import list2cmdline
import core.bundles
import core.modules.basic_modules
import core.modules.module_registry
import os
import os.path
import platform
import shutil
import sys

from PyQt4 import QtCore, QtGui
from packages.spreadsheet.basic_widgets import SpreadsheetCell, CellLocation
from packages.spreadsheet.spreadsheet_cell import QCellWidget
import packages.spreadsheet

from packages.controlflow.list_module import ListOfElements

import pyalps

basic = core.modules.basic_modules

config = ConfigurationObject()

##############################################################################
def _get_path(binary_file):
   
    if config.check('alpspath'):
        exename = os.path.join(config.alpspath, binary_file)
    else:
        exename = binary_file
    if platform.system()=='Windows':
        return exename +'.exe'
    else:
        return exename

def _get_tool_path(binary_file):
    if config.check('toolpath'):
        exename = os.path.join(config.toolpath, binary_file)
    else:
        exename = binary_file
    if platform.system()=='Windows':
        return exename +'.exe'
    else:
        return exename

def _get_default_mpi_procs():
    if config.check('mpiprocs'):
      return config.mpiprocs
    else:
      return 0

def _get_mpi_run():
    if config.check('mpirun'):
      cmd = "res =" + config.mpirun
      exec cmd
      return res
    else:
      return ['mpirun','-np']
      
class SystemCommand(Module):
    def execute(self,cmdline):
        if pyalps.executeCommand(cmdline) <> 0:
           raise ModuleError(self, 'Execution failed')

class SystemCommandLogged(Module):
    def execute(self,cmdline):
        logfile = self.interpreter.filePool.create_file(suffix='.log')
        result = pyalps.executeCommandLogged(cmdline,logfile.name)
        self.setResult('log_file', logfile)  
        if result <> 0:
          if platform.system()=='Darwin':
            pyalps.executeCommand(['open',logfile.name])
)         raise ModuleError(self, 'Execution failed')
    _output_ports = [('log_file',[basic.File])]


class OpenHTML(NotCacheable, SystemCommand):
    """ open the file using the system open command """
    def compute(self):
        cmdlist = 'false'
        if platform.system() == 'Windows':
          cmdlist = ['start','C:\Program Files\Internet Explorer\iexplore.exe']
        if platform.system()=='Darwin':
          cmdlist = ['open', '-a', 'Safari']
        if self.hasInputFromPort('file'):
           cmdlist += [self.getInputFromPort('file').name]
        if self.hasInputFromPort('files'):
           cmdlist += self.getInputFromPort('files')
        self.execute(cmdlist)
    _input_ports = [('file', [basic.File]),
                    ('files', [ListOfElements])]
    _input_ports = [('file', [basic.File]),
                    ('files', [ListOfElements])]

class TextFile(basic.Module):
    def compute(self):
        if self.hasInputFromPort('suffix'):
          f = self.interpreter.filePool.create_file(suffix=self.getInputFromPort('suffix'))
        else:
          f = self.interpreter.filePool.create_file()
        out = file(f.name,'w')
        out.write(self.getInputFromPort('text'))
        out.close()
        self.setResult('file',f)
    _input_ports = [('text', [basic.String]),
                    ('suffix', [basic.String])]
    _output_ports = [('file',[basic.File])]


class TextCell(SpreadsheetCell):
    """
    TextCell is a custom Module to view plain text files
    
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
  TextCellWidget has a QTextEdit to display plain text files
    
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



def dict_conf(l):
    if l[0] == '{' and l[-1] == '}':
        cmd = 'temp = ' + l
        exec cmd
        return temp
    else:
        pairs = l.split(',')
        temp = {}
        for pair in pairs:
            [k,v] = pair.split('=')
            temp[eval(k)] = eval(v)
        return temp

def dict_compute(self):
    if self.hasInputFromPort('value'):
        inps = self.forceGetInputListFromPort('value')
        result = {}
        for inp in inps:
            result.update(inp)
        self.setResult('value',result)
        self.setResult('value_as_string',str(result))


Dictionary = basic.new_constant('Dictionary', staticmethod(dict_conf), {},\
    staticmethod(lambda x: type(x) == dict))

def initialize(): pass


def selfRegister():

    reg = core.modules.module_registry.get_module_registry()

    reg.add_module(SystemCommand,namespace="Tools",abstract=True)
    reg.add_module(SystemCommandLogged,namespace="Tools",abstract=True)
    reg.add_module(OpenHTML,namespace="Tools")
    reg.add_module(TextFile,namespace="Tools")

    reg.add_module(TextCell,namespace="Tools")
    reg.add_input_port(TextCell, "Location", CellLocation)
    reg.add_input_port(TextCell, "File", basic.File)


    Dictionary.compute = dict_compute
    basic.init_constant(Dictionary)
#    reg.add_module(Dictionary,namespace="Parameters")
