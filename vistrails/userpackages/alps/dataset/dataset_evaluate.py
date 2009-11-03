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

			self.setResult('value',DataSets(d))

class EvalExpression_1to1(NotCacheable, Module):
	my_input_ports = [
		PortDescriptor("input",DataSets),
		PortDescriptor("source",basic.String,use_python_source=True)
	]
	my_output_ports = [
		PortDescriptor("output",DataSets)
	]

	def compute(self):
		if self.hasInputFromPort('input') and self.hasInputFromPort('source'):
			q = copy.deepcopy(self.getInputFromPort('input').sets)
			for s in q:
				x = s.x
				y = s.y

				code = self.getInputFromPort('source')
				proc_code = urllib.unquote(str(code))
				exec proc_code

				s.x = x
				s.y = y

			self.setResult('output',DataSets(q))

class EvalExpression_2to1(NotCacheable, Module):
	my_input_ports = [
		PortDescriptor("input1",DataSets),
		PortDescriptor("input2",DataSets),
		PortDescriptor("source",basic.String,use_python_source=True)
	]
	my_output_ports = [
		PortDescriptor("output",DataSets)
	]
	
	def compute(self):
		if self.hasInputFromPort('input1') and self.hasInputFromPort('input2') and self.hasInputFromPort('source'):
			q1 = copy.deepcopy(self.getInputFromPort('input1').sets)
			q2 = copy.deepcopy(self.getInputFromPort('input2').sets)
			
			results = []
			
			assert(len(q1) == len(q2))
			
			for iset in range(0,len(q1)):
				x1 = q1[iset].x
				y1 = q1[iset].y
				
				x2 = q2[iset].x
				y2 = q2[iset].y
				
				code = self.getInputFromPort('source')
				proc_code = urllib.unquote(str(code))
				exec proc_code
				
				result = DataSet()
				result.x = x
				result.y = y
				
				results.append(result)
			
			self.setResult('output',DataSets(results))

class EvalExpression_Allto1(NotCacheable, Module):
	my_input_ports = [
		PortDescriptor("input",DataSets),
		PortDescriptor("source",basic.String,use_python_source=True)
	]
	my_output_ports = [
		PortDescriptor("output",DataSets)
	]

	def compute(self):
		if self.hasInputFromPort('input') and self.hasInputFromPort('source'):
			x = np.array([])
			y = np.array([])
			
			res = DataSet()
			
			for s in self.getInputFromPort('input').sets:
				code = self.getInputFromPort('source')
				proc_code = urllib.unquote(str(code))
				exec proc_code
			
			self.setResult('output',DataSets(res))

class Select(Module):
	my_input_ports = [
		PortDescriptor("input",DataSets),
		PortDescriptor("source",basic.String,use_python_source=True)
	]
	my_output_ports = [
		PortDescriptor("output",DataSets)
	]

	def compute(self):
		if self.hasInputFromPort('input') and self.hasInputFromPort('source'):
			q = copy.deepcopy(self.getInputFromPort('input').sets)
			kept_sets = []
			for s in q:
				keep = True
				x = s.x
				y = s.y
				props = s.props

				code = self.getInputFromPort('source')
				proc_code = urllib.unquote(str(code))
				exec proc_code

				if keep:
					kept_sets.append(s)

			self.setResult('output',DataSets(kept_sets))
	
