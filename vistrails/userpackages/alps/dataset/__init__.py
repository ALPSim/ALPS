import core.modules.module_registry
import core.modules.basic_modules as basic
from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.modules.python_source_configure import PythonSourceConfigurationWidget

import urllib, copy
import numpy as np

# version = "0.0.3"
# name = "dataset"
# identifier = "org.comp-phys.alps.dataset"

# import sys, os
# try:
#     sys.path.append(os.environ['HOME'] + '/.vistrails/userpackages/alps/pyalps')
#     sys.path.append(os.environ['HOME'] + '/.vistrails/userpackages/alps/util')
# except KeyError:
#     raise EnvironmentError('Cannot find $HOME - do we live on Windows?')

from dataset_exceptions import *
from dataset_core import *
from dataset_evaluate import *
from dataset_fit import *
from dataset_load import *
from dataset_plot import *
from dataset_tools import *
from dataset_select import *

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
    
    if Descriptor in m.__bases__ or Selector in m.__bases__ or m == Selector:
        reg.add_output_port(m, 'output', (m,''))

def initialize():
    # We'll first create a local alias for the module registry so that
    # we can refer to it in a shorter way.
    reg = core.modules.module_registry.registry

    register(DataSets,'DataSet',abst=True)
    register(ConcatenateDataSets,'DataSet')
    register(ResultFiles,'DataSet',abst=True)
    
    register(ConstantDataSet,'DataSet|Load')
    register(GenerateDataSet,'DataSet|Load')
    register(LoadDataSet,'DataSet|Load')
    register(LoadAlpsFromTxt,'DataSet|Load')
    register(CustomLoader,'DataSet|Load')
    register(CollectXY,'DataSet|Load')
    register(LoadProperties,'DataSet|Load')
    register(LoadAlpsHdf5,'DataSet|Load')

    register(Transform,'DataSet|Evaluate')
    AddDataSetsInputPorts(TransformN, 5)
    register(TransformN,'DataSet|Evaluate')
    register(Reduce,'DataSet|Evaluate')
    register(GeneralTransform,'DataSet|Evaluate')
    
    register(AxisDescriptor,'DataSet|Plot')
    register(LegendDescriptor,'DataSet|Plot')
    register(Plotter,'DataSet|Plot',abst=True)
    register(PlotDescriptor,'DataSet|Plot')
    register(MplXYPlot,'DataSet|Plot')
    reg.add_module(PlotAsText,namespace='DataSet|Plot')
    reg.add_module(GraceXYPlot,namespace='DataSet|Plot')
    reg.add_module(GnuplotXYPlot,namespace='DataSet|Plot')
    
    register(FitPrototype,'DataSet|Fit',abst=True)
    register(PolyFit,'DataSet|Fit')
    register(NonlinearFit,'DataSet|Fit')
    
    register(SortByX,'DataSet|Tools')
    register(SelectXRange,'DataSet|Tools')
    register(WriteTxt,'DataSet|Tools')
    
    register(Selector,'DataSet|Select',abst=True)
    register(PropertySelector,'DataSet|Select')
    register(PropertyRangeSelector,'DataSet|Select')
    register(And,'DataSet|Select')
    register(Or,'DataSet|Select')
    register(Select,'DataSet',abst=True)
    register(Select,'DataSet|Select')
    register(SelectFiles,'DataSet|Select')
