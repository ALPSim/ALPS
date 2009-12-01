import core.modules.module_registry
import core.modules.basic_modules as basic
from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.modules.python_source_configure import PythonSourceConfigurationWidget

import urllib, copy
import numpy as np

from pyalps.util.dataset import DataSet

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

class ResultFiles(Module):
    my_input_ports = []
    my_output_ports = []
    
    def __init__(self):
        Module.__init__(self)

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
