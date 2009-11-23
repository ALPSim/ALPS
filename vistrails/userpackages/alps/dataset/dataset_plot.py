import core.modules.module_registry
import core.modules.basic_modules as basic
from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.modules.python_source_configure import PythonSourceConfigurationWidget
from packages.controlflow.list_module import ListOfElements

import urllib, copy
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

from dataset_exceptions import *
from dataset_core import *

class AxisDescriptor(Descriptor, Module):
	my_input_ports = [
		PortDescriptor('label',basic.String),
		PortDescriptor('min',basic.Float),
		PortDescriptor('max',basic.Float),
		PortDescriptor('logarithmic',basic.Boolean)
	]

class LegendDescriptor(Descriptor, Module):
	my_input_ports = [
		PortDescriptor('location',basic.Integer)
	]

class PlotDescriptor(Descriptor, Module):
	my_input_ports = [
		PortDescriptor('xaxis',AxisDescriptor),
		PortDescriptor('yaxis',AxisDescriptor),
		PortDescriptor('legend',LegendDescriptor),
		PortDescriptor('data',DataSets),
		PortDescriptor('title',basic.String)
	]

class MplXYPlot(Module):
	my_input_ports = [
		PortDescriptor('plot',PlotDescriptor),
		PortDescriptor('hide_buttons',basic.Boolean),
		PortDescriptor('source',basic.String)
	]
	my_output_ports = [PortDescriptor('unused',basic.String)]
	
	def __init__(self):
		Module.__init__(self)		
		self.colors = ['k','b','g','m','c','y']
	
	def hifp(self,m):
		return self.hasInputFromPort(m)
	def gifp(self,m):
		return self.getInputFromPort(m)
	
	def draw_lines(self):
		self.lines = []
		self.icolor = 0
		
		xlog = False
		ylog = False
		if 'xaxis' in self.plt and 'logarithmic' in self.plt['xaxis']:
			xlog = self.plt['xaxis']['logarithmic']
		if 'yaxis' in self.plt and 'logarithmic' in self.plt['yaxis']:
			ylog = self.plt['yaxis']['logarithmic']
		
		print self.plt	
		for q in self.plt['data']:
			if xlog and ylog:
				self.lines.append(plt.loglog(q.x,q.y,self.colors[self.icolor]))
			elif xlog:
				self.lines.append(plt.semilogx(q.x,q.y,self.colors[self.icolor]))
			elif ylog:
				self.lines.append(plt.semilogy(q.x,q.y,self.colors[self.icolor]))
			else:
				self.lines.append(plt.plot(q.x,q.y,self.colors[self.icolor]))
			
			if 'label' in q.props and q.props['label'] != 'none':
				self.lines[-1][0].set_label(q.props['label'])
			elif 'filename' in q.props:
				self.lines[-1][0].set_label(q.props['filename'])
			
			self.icolor = (self.icolor+1)%len(self.colors)
			
	def compute(self):
		if self.hifp('plot'):
			self.plt = self.gifp('plot')
		else:
			raise EmptyInputPort('plot')
		
		print self.plt
		
		self.draw_lines()
		
		if 'xaxis' in self.plt:
			if 'label' in self.plt['xaxis']:
				plt.xlabel(self.plt['xaxis']['label'])
			if 'min' in self.plt['xaxis'] and 'max' in self.plt['xaxis']:
				plt.xlim(self.plt['xaxis']['min'],self.plt['xaxis']['max'])
				
		if 'yaxis' in self.plt:
			if 'label' in self.plt['yaxis']:
				plt.ylabel(self.plt['yaxis']['label'])
			if 'min' in self.plt['yaxis'] and 'max' in self.plt['yaxis']:
				plt.ylim(self.plt['yaxis']['min'],self.plt['yaxis']['max'])
		
		if 'legend' in self.plt:
			if 'location' in self.plt['legend']:
				plt.legend(loc=self.plt['legend']['location'])
			else:
				plt.legend()
		
		if self.hifp('hide_buttons') and self.gifp('hide_buttons') == True:
			plt.get_current_fig_manager().toolbar.hide()
		
		if self.hasInputFromPort('source'):
			code = self.getInputFromPort('source')
			exec urllib.unquote(str(code))
		

class Plotter(NotCacheable, Module):
	my_input_ports = [
		PortDescriptor('data',DataSets),
		PortDescriptor('title',basic.String,hidden=True),
		PortDescriptor('xaxis',AxisDescriptor),
		PortDescriptor('yaxis',AxisDescriptor),
		PortDescriptor('legend',LegendDescriptor),
		PortDescriptor('hide_buttons',basic.Boolean,hidden=True),
		PortDescriptor('source',basic.String,use_python_source=True)
	]
	my_output_ports = [PortDescriptor('unused',basic.String)]

	def hifp(self,m):
		return self.hasInputFromPort(m)
	def gifp(self,m):
		return self.getInputFromPort(m)

	def compute(self):
		self.colors = ['k','b','g','m','c','y']
		
		if self.hifp('data'):
			s = self.gifp('data')
			lines = []
			icolor = 0
			for q in s:
				lines.append(plt.plot(q.x,q.y,self.colors[icolor]))
				icolor = (icolor+1)%len(self.colors)

				if q.props.has_key('label'):
					label = q.props['label']
					if label != 'none':
						lines[-1][0].set_label(q.props['label'])
				elif q.props.has_key('filename'):
					lines[-1][0].set_label(q.props['filename'])

			if self.hifp('title'):
				plt.title(self.gifp('title'))
			
			if self.hifp('xaxis'):
				xaxis = self.gifp('xaxis')
				if 'label' in xaxis:
					plt.xlabel(xaxis['label'])
				
				if 'min' in xaxis and 'max' in xaxis:
					if xaxis['min'] != xaxis['max']:
						plt.ylim(xaxis['min'],xaxis['max'])
			
			if self.hifp('yaxis'):
				yaxis = self.gifp('yaxis')
				if 'label' in yaxis:
					plt.ylabel(yaxis['label'])
				
				if 'min' in yaxis and 'max' in yaxis:
					if yaxis['min'] != yaxis['max']:
						plt.ylim(yaxis['min'],yaxis['max'])
			
			if self.hifp('legend'):
				legend = self.gifp('legend')
				if 'location' in legend:
					plt.legend(loc=legend['location'])
				else:
					plt.legend()
			
			if self.hifp('hide_buttons') and self.gifp('hide_buttons') == True:
				plt.get_current_fig_manager().toolbar.hide()
			
			if self.hasInputFromPort('source'):
				code = self.getInputFromPort('source')
				exec urllib.unquote(str(code))

			self.setResult('source','foo')
		else:
			raise EmptyInputPort('data')
