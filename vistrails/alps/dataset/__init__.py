# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 1994-2009 by Bela Bauer <bauerb@phys.ethz.ch>
# 
# This software is part of the ALPS libraries, published under the ALPS
# Library License; you can use, redistribute it and/or modify it under
# the terms of the license, either version 1 or (at your option) any later
# version.
#  
# You should have received a copy of the ALPS Library License along with
# the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
# available from http://alps.comp-phys.org/.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
# FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
# 
# ****************************************************************************

import core.modules.module_registry
import core.modules.basic_modules as basic
from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.modules.python_source_configure import PythonSourceConfigurationWidget

import urllib, copy
import numpy as np

# version = "0.0.3"
# name = "dataset"
# identifier = "org.comp-phys.alps.dataset"

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
        if input_port.porttype == SelftypePlaceholder:
            input_port.porttype = m
        
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
        if output_port.porttype == SelftypePlaceholder:
            output_port.porttype = m
        
        if output_port.hidden:
            reg.add_output_port(m, output_port.name, (output_port.porttype, output_port.description), True)
        else:
            reg.add_output_port(m, output_port.name, (output_port.porttype, output_port.description))
    
    if Descriptor in m.__bases__ or Selector in m.__bases__ or m == Selector:
        reg.add_output_port(m, 'output', (m,''))

def initialize(): pass

def selfRegister():

    # We'll first create a local alias for the module registry so that
    # we can refer to it in a shorter way.
    reg = core.modules.module_registry.registry

    register(DataSets,'DataSet',abst=True)
    register(ConcatenateDataSets,'DataSet')
    register(ResultFiles,'DataSet')
    
    register(ConstantDataSet,'DataSet|Load')
    register(GenerateDataSet,'DataSet|Load')
    register(LoadDataSet,'DataSet|Load')
    register(CustomLoader,'DataSet|Load')
    register(CollectXY,'DataSet|Load')
    register(LoadProperties,'DataSet|Load')
    register(LoadAlpsHdf5,'DataSet|Load')
    register(LoadSpectrumHdf5,'DataSet|Load')
    register(LoadBinningAnalysis,'DataSet|Load')
    register(LoadAlpsDiagData,'DataSet|Load')

    register(Transform,'DataSet|Evaluate')
    register(TransformProperties,'DataSet|Evaluate')
    AddDataSetsInputPorts(TransformN, 5)
    register(TransformN,'DataSet|Evaluate')
    register(Reduce,'DataSet|Evaluate')
    register(GeneralTransform,'DataSet|Evaluate')
    
    register(AxisDescriptor,'DataSet|Plot')
    register(LegendDescriptor,'DataSet|Plot')
    register(Plotter,'DataSet|Plot',abst=True)
    register(PlotDescriptor,'DataSet|Plot')
    register(MplXYPlot,'DataSet|Plot')
    reg.add_module(Convert2Text,namespace='DataSet|Plot')
    reg.add_module(GraceXYPlot,name='Convert2Grace',namespace='DataSet|Plot')
#    reg.add_module(GraceXYPlot,namespace='DataSet|Plot',abstract=True)
    reg.add_module(Convert2Gnuplot,namespace='DataSet|Plot')
    reg.add_module(LoadXMLPlot,namespace='DataSet|Plot')
    
    register(FitPrototype,'DataSet|Fit',abst=True)
    register(PolyFit,'DataSet|Fit')
    register(NonlinearFit,'DataSet|Fit')
    
    register(SortByX,'DataSet|Tools')
    register(SelectXRange,'DataSet|Tools')
    register(WriteTxt,'DataSet|Tools')
    register(CacheErasure,'DataSet|Tools')
    register(SetLabels,'DataSet|Tools')
    register(PrepareDictionary,'Tools')
    register(MakeScatter,'DataSet|Plot')
    
    register(Selector,'DataSet|Select',abst=True)
    register(PropertySelector,'DataSet|Select')
    register(PropertyRangeSelector,'DataSet|Select')
    register(ObservableSelector,'DataSet|Select')
    register(And,'DataSet|Select')
    register(Or,'DataSet|Select')
#    register(Select,'DataSet',abst=True)
    register(Select,'DataSet|Select')
    register(SelectFiles,'DataSet|Select')
    
    register(GroupBy,'DataSet|Hierarchy')
    register(GroupedTransform,'DataSet|Hierarchy')
    register(Flatten,'DataSet|Hierarchy')
    register(PrintHierarchyStructure,'DataSet|Hierarchy')
