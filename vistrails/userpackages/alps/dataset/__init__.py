import core.modules.module_registry
import core.modules.basic_modules as basic
from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.modules.python_source_configure import PythonSourceConfigurationWidget

import urllib, copy
import numpy as np

# version = "0.0.3"
# name = "dataset"
# identifier = "org.comp-phys.alps.dataset"

from dataset_core import *
from dataset_evaluate import *
from dataset_fit import *

def warn(m):
	print m

def register(m,ns,abst=False):
	reg = core.modules.module_registry.registry
	
	ups = False
	for input_port in m.my_input_ports:
		if input_port.use_python_source:
			ups = True
	if ups:
		reg.add_module(m, namespace=ns, configureWidgetType=PythonSourceConfigurationWidget, abstract=abst)
	else:
		reg.add_module(m, namespace=ns, abstract=abst)
	
	for input_port in m.my_input_ports:
		if ups and input_port.hidden:
			warn("use_python_source and hidden ports is currently not possible")
			if input_port.use_python_source:
				reg.add_input_port(m, input_port.name, (input_port.porttype, input_port.description), True)
			else:
				reg.add_input_port(m, input_port.name, (input_port.porttype, input_port.description))
		elif input_port.use_python_source or input_port.hidden:
			reg.add_input_port(m, input_port.name, (input_port.porttype, input_port.description), True)
		else:
			reg.add_input_port(m, input_port.name, (input_port.porttype, input_port.description))
	
	for output_port in m.my_output_ports:
		if output_port.hidden:
			reg.add_output_port(m, output_port.name, (output_port.porttype, output_port.description), True)
		else:
			reg.add_output_port(m, output_port.name, (output_port.porttype, output_port.description))

def initialize():
	# We'll first create a local alias for the module registry so that
	# we can refer to it in a shorter way.
	reg = core.modules.module_registry.registry

	basic.init_constant(Dictionary)
	register(Parameter,'DataSet')
	register(DataSets,'DataSet')
	register(ConcatenateDataSets,'DataSet')
	register(Select,'DataSet')
	
	register(ConstantDataSet,'DataSet::Load')
	register(LoadDataSet,'DataSet::Load')
	register(LoadAlpsFromTxt,'DataSet::Load')
	
	register(EvalExpression_1to1,'DataSet::Evaluate')
	register(EvalExpression_2to1,'DataSet::Evaluate')
	register(Plotter,'DataSet::Plot')
	
	register(FitPrototype,'DataSet::Fit',abst=True)
	register(PolyFit,'DataSet::Fit')
