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

class Loader:
	def __init__(self,filename,label,xcolumn,ycolumns,props={}):
		self.sets = []
		self.read_set(filename,label,xcolumn,ycolumns,props)
	
	def __init__(self):
		self.sets = []
	
	def read_set(self,filename,label,xcolumn,ycolumns,props={}):
		raw = np.loadtxt(filename).transpose()
		
		if len(raw.shape) == 0:
			res = DataSet()
			res.x = np.array([0])
			res.y = np.array([raw])
			res.props.update(props)
			res.props.update({'column':0,'label':copy.deepcopy(label),'filename':filename})
			self.sets.append(res)
			return
		
		for iyc in ycolumns:
			res = DataSet()
			res.props['filename'] = filename
			
			if xcolumn >= 0:
				res.x = raw[xcolumn]
				res.y = raw[iyc]
			else:
				if len(raw.strides) == 1:
					res.x = np.arange(0,len(raw))
					res.y = raw
				else:
					res.y = raw[iyc]
					res.x = np.arange(0,len(res.y))
			
			res.props['column'] = iyc
			res.props['label'] = copy.deepcopy(label)
			res.props.update(props)
			
			self.sets.append(res)

class LoadDataSet(Module):
	"""Load dataset from text data. Description of input ports:
	@file: The file, of which only the name is used.
	@label: A possible label for the dataset
	@x-column: The column of the data used for x data. Default: 0
	@y-columns: A ListOfElements of the columns. Default: 1
	Catch: if the file only has one columns, set x-column to -1!
	The following properties are set: filename, [label], column"""
	my_input_ports = [
		PortDescriptor("file",basic.File),
		PortDescriptor("label",basic.String),
		PortDescriptor("x-column",basic.Integer),
		PortDescriptor("y-columns",ListOfElements)
	]

	my_output_ports = [
		PortDescriptor("data",DataSets)
	]

	def compute(self):
		if self.hasInputFromPort('file'):
			f = self.getInputFromPort('file')
			filename = f.name
			print filename
			
			xc = 0
			yc = [1]
			if self.hasInputFromPort('x-column'):
				xc = int(self.getInputFromPort('x-column'))
				if xc < 0:
					yc = [0]
			if self.hasInputFromPort('y-columns'):
				yc = self.getInputFromPort('y-columns')
				print type(yc)
				print yc
				yc = [int(x) for x in yc]
			
			label = ''
			if self.hasInputFromPort('label'):
				label = self.getInputFromPort('label')
			
			l = Loader(filename,label,xc,yc)
			self.setResult('data', l.sets)

class CustomLoader(Module):
	my_input_ports = [
		PortDescriptor("source",basic.String,use_python_source=True),
		PortDescriptor("base_path",basic.String)
	]
	my_output_ports = [
		PortDescriptor("data",DataSets)
	]
	
	def compute(self):
		if self.hasInputFromPort('source'):
			code = self.getInputFromPort('source')
			proc_code = urllib.unquote(str(code))
			exec proc_code

class LoadAlpsFromTxt(Module):
	my_input_ports = [
		PortDescriptor('file',basic.File),
		PortDescriptor('parameter_name',basic.String)
	]
	my_output_ports = [
		PortDescriptor('data',DataSets)
	]

	def compute(self):
		if self.hasInputFromPort('file') and self.hasInputFromPort('parameter_name'):
			parname = self.getInputFromPort('parameter_name')
			filename = self.getInputFromPort('file').name
			lines = open(filename).readlines()
			par = 0
			sets = []

			data = []
			lines.append(parname+'=123') # evil hack
			for line in lines:
				if line.count(parname+'=') > 0:
					if len(data) > 0:
						res = DataSet()
						[res.x, res.y] = np.array(data).transpose()
						res.props[parname] = par
						res.props['label'] = parname + ' = ' + str(par)
						sets.append(copy.deepcopy(res))
					par = float(line.split('=')[1])
					data = []
				else:
					spl = line.split()
					idx = float(spl[0])
					val = float(spl[1])
					data.append([idx,val])

			self.setResult('data', sets)
		else:
			# throw something
			pass

class LoadAlpsHdf5(Module):
	my_input_ports = [
		PortDescriptor('file',basic.File),
		PortDescriptor('Measurements',ListOfElements)
	]
	
	my_output_ports = [
		PortDescriptor('data',DataSets)
	]	
	
	# Pre: file is a h5py file descriptor
	# Post: returns a list of all measurements saved in the file
	def FindAllMeasurements(self,file):
		return []
	
	# Pre: file is a h5py file descriptor
	# Post: returns DataSet with all parameters set
	def ReadMeasurementFromFile(self,file,measname):
		pass
	
	def compute(self):
		if self.hasInputFromPort('Measurements'):
			warn("Not implemented")
		else:
			# open file
			
			sets = []
			
			LOM = FindAllMeasurements(fd)
			for m in LOM:
				sets.append(ReadVariableFromFile(fd,m))
			
			self.setResult('data', sets)

class CollectXY(Module):
	my_input_ports = [
		PortDescriptor('for-each',ListOfElements),
		PortDescriptor('observable',basic.String),
		PortDescriptor('input',DataSets)
	]
	my_output_ports = [
		PortDescriptor('output',DataSets)
	]
	
	def compute(self):
		if self.hasInputFromPort('for-each') and self.hasInputFromPort('observable'):
			# find all possible values for each for-each
			sets = self.getInputFromPort('input')
			for_each = self.getInputFromPort('for-each')
			observable = self.getInputFromPort('observable')
			
			for_each_sets = {}
			for iset in sets:
				if iset.props['observable'] != observable:
					continue
				
				fe_par_set = []
				for m in for_each:
					fe_par_set.append(iset.props[m])
				
				if fe_par_set in for_each_sets:
					for_each_sets[fe_par_set].append(iset)
				else:
					for_each_sets[fe_par_set] = [iset]
			
			for k,v in for_each_sets.items():
				res = DataSet()
				res.props = v[0].props
				for im in range(0,for_each):
					m = for_each[im]
					res.props[m] = k[im]
				
				for x in v:
					res.x = np.concatenate(res.x, x.x)
					res.y = np.concatenate(res.y, x.y)
				
				order = np.argsort(res.x)
				res.x = res.x[order]
				res.y = res.y[order]
				
				for_each_sets[k] = res
			
			self.setResult('output',for_each_sets.values())
				
		else:
			raise EmptyInputPort('for-each || observable')
