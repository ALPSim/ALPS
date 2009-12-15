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
from dataset_fit import *

class SortByX(FitPrototype):
	def transform(self,data):
		order = np.argsort(data.x)
		data.x = data.x[order]
		data.y = data.y[order]

class SelectXRange(FitPrototype):
    my_input_ports = [
        PortDescriptor('min',basic.Float),
        PortDescriptor('max',basic.Float)
    ]
    
    def transform(self,data):
        min = self.getInputFromPort('min')
        max = self.getInputFromPort('max')
        
        selection = (data.x >= min) & (data.x <= max)
        data.x = data.x[selection]
        data.y = data.y[selection]

class WriteTxt(Module):
	my_input_ports = [PortDescriptor('input',DataSets)]
	my_output_ports = []
	
	def compute(self):
		if self.hasInputFromPort('input'):
			for s in self.getInputFromPort('input'):
				if 'filename' in s.props:
					data = np.array([s.x,s.y]).transpose()
					np.savetxt(s.props['filename'],data)

class CacheErasure(NotCacheable,Module):
    my_input_ports = [PortDescriptor('input',DataSets)]
    my_output_ports = [PortDescriptor('output',DataSets)]
    
    def compute(self):
        self.setResult('output', self.getInputFromPort('input'))
    
