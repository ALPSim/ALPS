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
				if 'label' in xaxis.data:
					print 'xlabel = ' + xaxis.data['label']
					plt.xlabel(xaxis.data['label'])
				
				if 'min' in xaxis.data and 'max' in xaxis.data:
					if xaxis.data['min'] != xaxis.data['max']:
						plt.ylim(xaxis.data['min'],xaxis.data['max'])
			
			if self.hifp('yaxis'):
				yaxis = self.gifp('yaxis')
				if 'label' in yaxis.data:
					print 'ylabel = ' + yaxis.data['label']
					plt.ylabel(yaxis.data['label'])
				
				if 'min' in yaxis.data and 'max' in yaxis.data:
					if yaxis.data['min'] != yaxis.data['max']:
						plt.ylim(yaxis.data['min'],yaxis.data['max'])
			
			if self.hifp('legend'):
				legend = self.gifp('legend')
				if 'location' in legend.data:
					plt.legend(loc=legend.data['location'])
			
			if self.hifp('hide_buttons') and self.gifp('hide_buttons') == True:
				plt.get_current_fig_manager().toolbar.hide()
			
			if self.hasInputFromPort('source'):
				code = self.getInputFromPort('source')
				exec urllib.unquote(str(code))

			self.setResult('source','foo')
		else:
			raise EmptyInputPort('data')
