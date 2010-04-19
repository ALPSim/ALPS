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
import copy

from dataset_core import *
from dataset_exceptions import *
from dataset_fit import *

from pyalps.hlist import deep_flatten, flatten, happly, hmap, depth
from pyalps.dict_intersect import dict_difference, dict_intersect

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
    
class SetLabels(Module):
    my_input_ports = [
        PortDescriptor('input',DataSets),
        PortDescriptor('label_props',ListOfElements)
    ]
    my_output_ports = [
        PortDescriptor('output',DataSets)
    ]
    
    # overwrite default behaviour to deepcopy everything
    def compute(self):
        q = self.getInputFromPort('input')
        labels = self.getInputFromPort('label_props')
        
        def f(x):
            ret = DataSet()
            ret.x = x.x
            ret.y = x.y
            ret.props = copy.deepcopy(x.props)
            labelstr = ''
            for label in labels:
                if label != labels[0]:
                    labelstr += ', '
                if type(x.props[label]) == str:
                    labelstr += '%s = %s' % (label,x.props[label])
                else:
                    labelstr += '%s = %.4s' % (label,x.props[label])
            ret.props['label'] = labelstr
            return ret
        
        q2 = hmap(f, q)
        self.setResult('output', q2)

class CycleColors(Module):
    my_input_ports = [
        PortDescriptor('input',DataSets),
        PortDescriptor('for-each',ListOfElements),
        PortDescriptor('colors',ListOfElements)
    ]
    my_output_ports = [
        PortDescriptor('output',DataSets)
    ]
    
    def compute(self):
        colors = ['k','b','g','m','c','y']
        if self.hasInputFromPort('colors'):
            colors = self.getInputFromPort('colors')
        
        input = self.getInputFromPort('input')
        foreach = self.getInputFromPort('for-each')
        
        all = {}
        for q in flatten(input):
            key = tuple([q.props[k] for k in foreach])
            all[key] = ''
        
        icolor = 0
        for k in all.keys():
            all[k] = colors[icolor]
            icolor = (icolor+1)%len(colors)
        
        for q in flatten(input):
            key = tuple([q.props[k] for k in foreach])
            q.props['color'] = all[key]
        
        self.setResult('output', input)

class CycleMarkers(Module):
    my_input_ports = [
        PortDescriptor('input',DataSets),
        PortDescriptor('for-each',ListOfElements),
        PortDescriptor('markers',ListOfElements)
    ]
    my_output_ports = [
        PortDescriptor('output',DataSets)
    ]
    
    def compute(self):
        markers = ['s', 'o', '^', '>', 'v', '<', 'd', 'p', 'h', '8', '+', 'x']
        if self.hasInputFromPort('markers'):
            markers = self.getInputFromPort('markers')
        
        input = self.getInputFromPort('input')
        foreach = self.getInputFromPort('for-each')
        
        all = {}
        for q in flatten(input):
            key = tuple([q.props[k] for k in foreach])
            all[key] = ''
        
        imarker = 0
        for k in all.keys():
            all[k] = markers[imarker]
            imarker = (imarker+1)%len(markers)
        
        for q in flatten(input):
            key = tuple([q.props[k] for k in foreach])
            q.props['marker'] = all[key]
        
        self.setResult('output', input)

class MakeScatter(FitPrototype):
    def transform(self,data):
        data.props['line'] = 'scatter'

class Flatten(Module):
    my_input_ports = [PortDescriptor('input',DataSets)]
    my_output_ports = [PortDescriptor('output',DataSets)]
    
    def compute(self):
        self.setResult('output', deep_flatten(self.getInputFromPort('input')))

class PrepareDictionary(Module):
    my_input_ports = [PortDescriptor('source',basic.String,use_python_source=True)]
    my_output_ports = [PortDescriptor('output',basic.Dictionary)]
    
    def compute(self):
        lines = self.getInputFromPort('source')
        lines = urllib.unquote(str(lines)).split('\n')
        lines1 = []
        for line in lines:
            if line.startswith('#'):
                continue
            else:
                lines1.append(line.strip())
        
        d = {}
        print lines1
        for line in lines1:
            pair = line.split('=')
            if len(pair) == 2:
                d[pair[0]] = pair[1]
        
        self.setResult('output', d)

class PrintHierarchyStructure(Module):
    my_input_ports = [PortDescriptor('input',DataSets)]
    my_output_ports = []
    
    def compute(self):
        q = self.getInputFromPort('input')
        
        all_props = [v.props for v in flatten(q)]
        diff_props = dict_difference(all_props)
        
        def f(x):
            d = {}
            for k in diff_props:
                d[k] = x.props[k]
            return d
        
        q2 = hmap(f, q)
        print q2
