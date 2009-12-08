import core.modules.module_registry
import core.modules.basic_modules as basic
from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.modules.python_source_configure import PythonSourceConfigurationWidget

from packages.controlflow.list_module import ListOfElements

import urllib, copy
import numpy as np

from pyalps.util.dataset import DataSet,ResultFile

class PortDescriptor:
    def __init__(self, name, porttype, description='', use_python_source=False, hidden=False):
        self.name = name
        self.porttype = porttype
        self.description = description
        self.use_python_source = use_python_source
        self.hidden = hidden
    
    name = ''
    porttype = basic.String
    description = ''
    use_python_source = False
    hidden = False

class DataSets(Module):
    my_input_ports = []
    my_output_ports = []

    def __init__(self):
        Module.__init__(self)

class SelftypePlaceholder:
    is_placeholder = True

class ResultFiles(Module):
    my_input_ports = [
        PortDescriptor('filenames', ListOfElements),
        PortDescriptor('resultfiles', SelftypePlaceholder)
    ]
    my_output_ports = [
        PortDescriptor('filenames', ListOfElements),
        PortDescriptor('resultfiles', SelftypePlaceholder)
    ]
    
    def __init__(self):
        Module.__init__(self)
    
    def compute(self):
        rf = []
        if self.hasInputFromPort('resultfiles'):
            rf += self.getInputFromPort('resultfiles')
        if self.hasInputFromPort('filenames'):
            rf += [ResultFile(x) for x in self.getInputFromPort('filenames')]
        
        self.setResult('resultfiles', rf)
        self.setResult('filenames',[x.props['filename'] for x in rf])

class ConcatenateDataSets(Module):
    """Takes the sets from several DataSets from the input port and puts them
    into one DataSets object"""
    my_input_ports = [
        PortDescriptor('input',DataSets)
    ]

    my_output_ports = [
        PortDescriptor('value',DataSets)
    ]

    def compute(self):
        if self.hasInputFromPort('input'):
            inps = self.forceGetInputListFromPort('input')
            sets = []
            for inp in inps:
                sets = sets + inp
            self.setResult('value',sets)

class Descriptor:
    my_input_ports = []
    my_output_ports = []
    
    def __init__(self):
        Module.__init__(self)
    
    def compute(self):
        ret = {}
        for ip in self.my_input_ports:
            if self.hasInputFromPort(ip.name):
                ret[ip.name] = self.getInputFromPort(ip.name)
        self.setResult('output',ret)
