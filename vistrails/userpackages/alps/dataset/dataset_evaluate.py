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

class LoadDataSet(Module):
	"""Load dataset from text data. Description of input ports:
	@file: The file, of which only the name is used.
	@label: A possible label for the dataset
	@x-column: The column of the data used for x data. Default: 0
	@y-columns: A ListOfElements of the columns. Default: 1
	Catch: if the file only has one columns, set x-column to -1!
	The following properties are set: filename, [labe], column"""
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
			filename = self.getInputFromPort('file').get_name()
			
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
			
			raw = np.loadtxt(filename).transpose()
			
			sets = []
			for iyc in yc:
				res = DataSet()
				res.props['filename'] = filename
				
				if xc >= 0:
					res.x = raw[xc]
					res.y = raw[iyc]
				else:
					if len(raw.strides) == 1:
						raw.x = np.arange(0,len(raw))
						raw.y = raw
					else:
						res.y = raw[yc]
						res.x = np.arange(0,len(res.y))
				
				res.props['column'] = iyc
				
				if self.hasInputFromPort('label'):
					res.props['label'] = copy.deepcopy(self.getInputFromPort('label'))
				sets.append(res)

			self.setResult('data', DataSets(sets))

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
			filename = self.getInputFromPort('file').get_name()
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

			self.setResult('data', DataSets(sets))

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

			self.setResult('output',DataSets(q))
	
class Plotter(NotCacheable, Module):
	my_input_ports = [
		PortDescriptor('data',DataSets),
		PortDescriptor('title',basic.String,hidden=True),
		PortDescriptor('xlabel',basic.String,hidden=True),
		PortDescriptor('ylabel',basic.String,hidden=True),
		PortDescriptor('x0',basic.Float,hidden=True),
		PortDescriptor('x1',basic.Float,hidden=True),
		PortDescriptor('y0',basic.Float,hidden=True),
		PortDescriptor('y1',basic.Float,hidden=True),
		PortDescriptor('display_legend',basic.Boolean,hidden=True),
		PortDescriptor('hide_buttons',basic.Boolean,hidden=True),
		PortDescriptor('source',basic.String,use_python_source=True)
	]
	my_output_ports = [PortDescriptor('unused',basic.String)]

	colors = ['k','b','g','m','c','y']

	def hifp(self,m):
		return self.hasInputFromPort(m)
	def gifp(self,m):
		return self.getInputFromPort(m)

	def compute(self):
		if self.hifp('data'):
			s = self.gifp('data')
			lines = []
			icolor = 0
			for q in s.sets:
				lines.append(plt.plot(q.x,q.y,self.colors[icolor]))
				icolor = (icolor+1)%len(self.colors)

				if q.props.has_key('label'):
					lines[-1][0].set_label(q.props['label'])
				elif q.props.has_key('filename'):
					lines[-1][0].set_label(q.props['filename'])

			if self.hifp('title'):
				plt.title(self.gifp('title'))
			if self.hifp('xlabel'):
				plt.xlabel(self.gifp('xlabel'))
			if self.hifp('ylabel'):
				plt.ylabel(self.gifp('ylabel'))
			if self.hifp('x0') and self.hifp('x1'):
				plt.xlim(self.gifp('x0'), self.gifp('x1'))
			if self.hifp('y0') and self.hifp('y1'):
				plt.ylim(self.gifp('y0'), self.gifp('y1'))
			if self.hifp('display_legend') and self.gifp('display_legend') == True:
				plt.legend()
			if self.hifp('hide_buttons') and self.gifp('hide_buttons') == True:
				plt.get_current_fig_manager().toolbar.hide()
			
			if self.hasInputFromPort('source'):
				code = self.getInputFromPort('source')
				exec urllib.unquote(str(code))

			self.setResult('source','foo')