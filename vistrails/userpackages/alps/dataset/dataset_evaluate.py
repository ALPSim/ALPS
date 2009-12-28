import core.modules.module_registry
import core.modules.basic_modules as basic
from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.modules.python_source_configure import PythonSourceConfigurationWidget
from packages.controlflow.list_module import ListOfElements

import urllib, copy
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

from dataset_core import *
from dataset_exceptions import *

class ConstantDataSet(Module):
    """Create a constant dataset and store into DataSets"""
    my_input_ports = [
        PortDescriptor("value",basic.Float),
        PortDescriptor("length",basic.Integer)
    ]

    my_output_ports = [
        PortDescriptor("value",DataSets)
    ]

    def compute(self):
        if self.hasInputFromPort('value') and self.hasInputFromPort('length'):
            value = self.getInputFromPort('value')
            length = self.getInputFromPort('length')

            d = DataSet()
            d.x = np.arange(0,length)
            d.y = value + 0*d.x

            self.setResult('value',[d])
        else:
            raise EmptyInputPort('value || source')

class GenerateDataSet(Module):
    my_input_ports = [PortDescriptor('source',basic.String,use_python_source=True)]
    my_output_ports = [PortDescriptor('output',DataSets)]
    
    def compute(self):
        if self.hasInputFromPort('source'):
            code = self.getInputFromPort('source')
            proc_code = urllib.unquote(str(code))
            exec proc_code
            
            try:
                if type(result) == list:
                    self.setResult('output',result)
                else:
                    self.setResult('output',[result])
            except NameError:
                raise InvalidInput("Generate result!")
                
        else:
            raise EmptyInputPort('source')

class Transform(Module):
    my_input_ports = [
        PortDescriptor("input",DataSets),
        PortDescriptor("source",basic.String,use_python_source=True)
    ]
    my_output_ports = [
        PortDescriptor("output",DataSets)
    ]
    deepcopy_all = True

    def compute(self):
        if self.hasInputFromPort('input') and self.hasInputFromPort('source'):
            if self.deepcopy_all:
                q = copy.deepcopy(self.getInputFromPort('input'))
            else:
                q = self.getInputFromPort('input')
                newsets = []
            for s in q:
                x = s.x
                y = s.y
                if self.deepcopy_all:
                    props = s.props
                else:
                    props = copy.deepcopy(s.props)

                code = self.getInputFromPort('source')
                proc_code = urllib.unquote(str(code))
                exec proc_code
                
                if self.deepcopy_all:
                    s.x = x
                    s.y = y
                else:
                    newsets.append(DataSet())
                    newsets[-1].x = x
                    newsets[-1].y = y
                    newsets[-1].props = props

            if self.deepcopy_all:
                self.setResult('output',q)
            else:
                self.setResult('output',newsets)
        else:
            raise EmptyInputPort('input || source')

class TransformProperties(Transform):
    deepcopy_all = False

def AddDataSetsInputPorts(m, Nmax):
    for i in range(0,Nmax):
        m.my_input_ports.append(PortDescriptor('input'+str(i),DataSets))

class TransformN(Module):
    my_input_ports = [
        PortDescriptor("source",basic.String,use_python_source=True)
    ]
    my_output_ports = [
        PortDescriptor("output",DataSets)
    ]
    
    def compute(self):
        if self.hasInputFromPort('source'):
            Nports = len(self.my_input_ports)-1
            print Nports
            
            inputs = []
            for i in range(0,Nports):
                port = 'input'+str(i)
                if self.hasInputFromPort(port):
                    r = self.getInputFromPort(port)
                    inputs.append(copy.deepcopy(r))
            Ninputs = len(inputs)
            
            results = []
            
            for i in range(1,Ninputs):
                if len(inputs[0]) != len(inputs[i]):
                    raise InvalidInput("Input lengths don't match")
            
            for iset in range(0,len(inputs[0])):
                for iport in range(0,Ninputs):
                    cmd = 'set' + str(iport) + ' = inputs[' + str(iport) + '][' + str(iset) + ']'
                    exec cmd
                
                result = DataSet()
                
                code = self.getInputFromPort('source')
                proc_code = urllib.unquote(str(code))
                exec proc_code
                
                results.append(result)
            
            self.setResult('output',results)
        else:
            raise EmptyInputPort('source')

class Reduce(Module):
    my_input_ports = [
        PortDescriptor("input",DataSets),
        PortDescriptor("source",basic.String,use_python_source=True)
    ]
    my_output_ports = [
        PortDescriptor("output",DataSets)
    ]

    def compute(self):
        if self.hasInputFromPort('input') and self.hasInputFromPort('source'):
            result = DataSet()
            
            for s in self.getInputFromPort('input'):
                code = self.getInputFromPort('source')
                proc_code = urllib.unquote(str(code))
                exec proc_code
            
            self.setResult('output',[result])
        else:
            raise EmptyInputPort('input || source')

class GeneralTransform(Module):
    my_input_ports = [
        PortDescriptor("input",DataSets),
        PortDescriptor("source",basic.String,use_python_source=True)
    ]
    my_output_ports = [
        PortDescriptor("output",DataSets)
    ]
    
    def compute(self):
        if self.hasInputFromPort('input') and self.hasInputFromPort('source'):
            code = self.getInputFromPort('source')
            proc_code = urllib.unquote(str(code))
            exec proc_code
        else:
            raise EmptyInputPort('input || source')
