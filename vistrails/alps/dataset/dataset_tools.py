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
from packages.controlflow.list_module import ListOfElements

import urllib, copy
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

from dataset_core import *
from dataset_exceptions import *
from dataset_fit import *

class SortByX(FitPrototype):
    def transform(self,data):
        order = np.argsort(data.x)
        data.x = data.x[order]
        data.y = data.y[order]

class SelectXRange(FitPrototype):
    my_input_ports = [
        PortDescriptor('min',basic.Float),
        PortDescriptor('max',basic.Float)
    ]
    
    def transform(self,data):
        min = self.getInputFromPort('min')
        max = self.getInputFromPort('max')
        
        selection = (data.x >= min) & (data.x <= max)
        data.x = data.x[selection]
        data.y = data.y[selection]

class WriteTxt(Module):
    my_input_ports = [PortDescriptor('input',DataSets)]
    my_output_ports = []
    
    def compute(self):
        if self.hasInputFromPort('input'):
            for s in self.getInputFromPort('input'):
                if 'filename' in s.props:
                    data = np.array([s.x,s.y]).transpose()
                    np.savetxt(s.props['filename'],data)

class CacheErasure(NotCacheable,Module):
    my_input_ports = [PortDescriptor('input',DataSets)]
    my_output_ports = [PortDescriptor('output',DataSets)]
    
    def compute(self):
        self.setResult('output', self.getInputFromPort('input'))
    
class SetLabels(FitPrototype):
    my_input_ports = [
        PortDescriptor('label_props',ListOfElements)
    ]
    
    # overwrite default behaviour to deepcopy everything
    def compute(self):
        FitPrototype.property_compute(self)
    
    def transform(self,data):
        labelstr = ''
        labels = self.getInputFromPort('label_props')
        for label in labels:
            if label != labels[0]:
                labelstr += ', '
            labelstr += '%s = %.4s' % (label,data.props[label])
        data.props['label'] = labelstr

class MakeScatter(FitPrototype):
    def transform(self,data):
        data.props['line'] = 'scatter'
