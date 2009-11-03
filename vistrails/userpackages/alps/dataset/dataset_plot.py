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

class AxisDescriptor(Descriptor, Module):
	my_input_ports = [
		PortDescriptor('label',basic.String),
		PortDescriptor('min',basic.Float),
		PortDescriptor('max',basic.Float)
	]

class LegendDescriptor(Descriptor, Module):
	my_input_ports = [
		PortDescriptor('location',basic.Integer)
	]

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
			
			if self.hifp('xaxis'):
				xaxis = self.gifp('xaxis')
				try:
					plt.xlabel(xaxis.data['label'])
				except KeyError:
					pass
				try:
					if xaxis.data['min'] != xaxis.data['max']:
						plt.xlim(xaxis.data['min'],xaxis.data['max'])
				except KeyError:
					pass
			
			if self.hifp('yaxis'):
				yaxis = self.gifp('yaxis')
				try:
					plt.ylabel(yaxis.data['label'])
				except KeyError:
					pass
				try:
					if yaxis.data['min'] != yaxis.data['max']:
						plt.ylim(yaxis.data['min'],yaxis.data['may'])
				except KeyError:
					pass
			
			if self.hifp('legend'):
				legend = self.gifp('legend')
				try:
					plt.legend(loc=legend.data['location'])
				except KeyError:
					plt.legend()
			
			if self.hifp('hide_buttons') and self.gifp('hide_buttons') == True:
				plt.get_current_fig_manager().toolbar.hide()
			
			if self.hasInputFromPort('source'):
				code = self.getInputFromPort('source')
				exec urllib.unquote(str(code))

			self.setResult('source','foo')