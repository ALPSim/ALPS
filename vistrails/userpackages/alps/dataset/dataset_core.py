import core.modules.module_registry
import core.modules.basic_modules as basic
from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.modules.python_source_configure import PythonSourceConfigurationWidget

import urllib, copy
import numpy as np

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

Dictionary.compute = dict_compute

class Parameter(Module):
	my_input_ports = [
		PortDescriptor('name',basic.String),
		PortDescriptor('value',basic.String)
	]
	
	my_output_ports = [
		PortDescriptor('value', Dictionary)
	]
	
	def compute(self):
		if self.hasInputFromPort('name') and self.hasInputFromPort('value'):
			name = eval(self.getInputFromPort('name'))
			value = eval(self.getInputFromPort('value'))
			self.setResult('value',{name: value})

class DataSet:
	def __deepcopy__(self,memo):
		ret = DataSet()
		ret.props = copy.deepcopy(self.props,memo)
		ret.x = copy.deepcopy(self.x,memo)
		ret.y = copy.deepcopy(self.y,memo)
		ret.dx = copy.deepcopy(self.dx,memo)
		ret.dy = copy.deepcopy(self.dy,memo)
		return ret

	def __init__(self):
		self.props = {}
		self.x = np.array([])
		self.y = np.array([])

		self.dx = np.array([])
		self.dy = np.array([])

class DataSets(Module):
	my_input_ports = []
	my_output_ports = []

	def __init__(self,sets_ = None):
		Module.__init__(self)
		self.sets = []
		if sets_ != None:
			if type(sets_) == type(self.sets):
				self.sets = sets_
			else:
				self.sets = [sets_]

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
				sets = sets + inp.sets
			self.setResult('value',DataSets(sets))

class Descriptor:
	my_input_ports = []
	my_output_ports = []
	
	def compute(self):
		for ip in self.my_input_ports:
			if self.hasInputFromPort(ip.name):
				self.data[ip.name] = self.getInputFromPort(ip.name)
		self.setResult('output',self)
	
	def __init__(self):
		Module.__init__(self)
		self.data = {}
