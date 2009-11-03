import core.modules.module_registry
import core.modules.basic_modules as basic
from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.modules.python_source_configure import PythonSourceConfigurationWidget

import urllib, copy
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize, polyfit
import fit_wrapper as fw

from dataset_core import *

class FitPrototype(Module):
	my_input_ports = [
		PortDescriptor("input",DataSets)
	]
	my_output_ports = [
		PortDescriptor("output",DataSets)
	]

	def compute(self):
		if self.hasInputFromPort('input'):
			q = copy.deepcopy(self.getInputFromPort('input').sets)
			for s in q:
				s = self.transform(s)

			self.setResult('output',DataSets(q))

class PolyFit(FitPrototype):
	my_input_ports = FitPrototype.my_input_ports + [PortDescriptor("degree",basic.Integer)]
	my_output_ports = FitPrototype.my_output_ports
	
	def transform(self, data):
		degree = 0
		if self.hasInputFromPort('degree'):
			degree = self.getInputFromPort('degree')
		else:
			degree = 1
		
		fit_parms = polyfit(data.x, data.y, degree)
		data.props['fit_parameters'] = fit_parms
		
		data.y = 0*data.x
		for deg in range(0,degree+1):
			data.y = data.y + fit_parms[deg]*data.x**(degree-deg)
		
		return data
	