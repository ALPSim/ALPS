import core.modules.module_registry
import core.modules.basic_modules as basic
from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.modules.python_source_configure import PythonSourceConfigurationWidget
from packages.controlflow.list_module import ListOfElements

import urllib, copy
import numpy as np

from dataset_core import *
from dataset_exceptions import *

class Selector(Module):
	my_input_ports = []
	my_output_ports = []
	
	def compute(self):
		for inp in self.my_input_ports:
			if not self.hasInputFromPort(inp.name):
				raise EmptyInputPort(inp.name)
		self.setResult('output',self)

class PropertySelector(Selector):
	my_input_ports = [
		PortDescriptor('property_name',basic.String),
		PortDescriptor('property_value',basic.String)
	]
	
	def decide(self,ds):
		pn = self.getInputFromPort('property_name')
		if pn in ds.props:
			pv = type(ds.props[pn])(self.getInputFromPort('property_value'))
			if ds.props[pn] == pv:
				return True
		return False

class PropertyRangeSelector(Selector):
	my_input_ports = [
		PortDescriptor('property_name',basic.String),
		PortDescriptor('min',basic.String),
		PortDescriptor('max',basic.String)
	]
	
	def decide(self,ds):
		pn = self.getInputFromPort('property_name')
		if pn in ds.props:
			pmin = type(ds.props[pn])(self.getInputFromPort('min'))
			pmax = type(ds.props[pn])(self.getInputFromPort('max'))
			if ds.props[pn] >= pmin and ds.props[pn] <= pmax:
				return True
		return False

class Select(Module):
	my_input_ports = [
		PortDescriptor("input",DataSets),
		PortDescriptor("source",basic.String,use_python_source=True),
		PortDescriptor('select',Selector)
	]
	my_output_ports = [
		PortDescriptor("kept",DataSets),
		PortDescriptor("discarded",DataSets)
	]

	def compute(self):
		if self.hasInputFromPort('input') and self.hasInputFromPort('source'):
			q = copy.deepcopy(self.getInputFromPort('input'))
			kept_sets = []
			disc_sets = []
			
			code = self.getInputFromPort('source')
			proc_code = urllib.unquote(str(code))
			
			cmd = 'def fn(x,y,props):\n'
			for line in proc_code.split('\n'):
				cmd = cmd + '\t' + line + '\n'
			exec cmd
				
			for s in q:
				if fn(s.x,s.y,s.props):
					kept_sets.append(s)
				else:
					disc_sets.append(s)
					
			self.setResult('kept',kept_sets)
			self.setResult('discarded',disc_sets)
		elif self.hasInputFromPort('input') and self.hasInputFromPort('select'):
			s = self.getInputFromPort('select')
			q = copy.deepcopy(self.getInputFromPort('input'))
			kept_sets = []
			disc_sets = []
			
			for iq in q:
				if s.decide(iq):
					kept_sets.append(iq)
				else:
					disc_sets.append(iq)
			
			self.setResult('kept',kept_sets)
			self.setResult('discarded',disc_sets)
		else:
			raise EmptyInputPort('input || source')

class And(Selector):
	my_input_ports = [
		PortDescriptor('selectors',Selector)
	]
	
	def compute(self):
		self.selectors = []
		if self.hasInputFromPort('selectors'):
			self.selectors = self.forceGetInputListFromPort('selectors')
		self.setResult('output',self)
		
	def decide(self,ds):
		for selector in self.selectors:
			if selector.decide(ds) == False:
				return False
		return True

class Or(And):
	def decide(self,ds):
		for selector in self.selectors:
			if selector.decide(ds):
				return True
		return False
