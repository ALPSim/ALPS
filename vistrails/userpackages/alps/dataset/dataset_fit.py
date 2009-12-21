import core.modules.module_registry
import core.modules.basic_modules as basic
from core.modules.vistrails_module import Module, ModuleError, NotCacheable
from core.modules.python_source_configure import PythonSourceConfigurationWidget
from packages.controlflow.list_module import ListOfElements

import urllib, copy
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize, polyfit
import pyalps.util.fit_wrapper as fw

from dataset_core import *
from dataset_exceptions import *

class FitPrototype(Module):
    my_input_ports = [
        PortDescriptor("input",DataSets)
    ]
    my_output_ports = [
        PortDescriptor("output",DataSets)
    ]

    def compute(self):
        if self.hasInputFromPort('input'):
            q = copy.deepcopy(self.getInputFromPort('input'))
            for s in q:
                s = self.transform(s)

            self.setResult('output',q)

class PolyFit(FitPrototype):
    my_input_ports = FitPrototype.my_input_ports + [PortDescriptor("degree",basic.Integer)]
    my_output_ports = FitPrototype.my_output_ports
    
    def transform(self, data):
        degree = 0
        if self.hasInputFromPort('degree'):
            degree = self.getInputFromPort('degree')
        else:
            degree = 1
        
        fit_parms = polyfit(data.x, data.y, degree)
        data.props['fit_parameters'] = fit_parms
        
        data.y = 0*data.x
        for deg in range(0,degree+1):
            data.y = data.y + fit_parms[deg]*data.x**(degree-deg)
        
        return data

class NonlinearFit(FitPrototype):
    my_input_ports = FitPrototype.my_input_ports + \
    [
        PortDescriptor('parameters',ListOfElements),
        PortDescriptor('source',basic.String,use_python_source=True),
        PortDescriptor('xrange_min',basic.Float),
        PortDescriptor('xrange_max',basic.Float)
    ]
    
    def transform(self, data):
        pars = []
        raw_pars = self.getInputFromPort('parameters')
        cmd = ''
        for i in range(0,len(raw_pars)):
            p = raw_pars[i][0]
            if p == 'p':
                raise InvalidInput('p is not a good parameter name')
            cmd += 'global ' + p + '\n'
            cmd += p + ' = fw.Parameter(raw_pars[' + str(i) + '][1])\n'
            cmd += 'pars.append(' + p + ')\n'
        exec cmd
        
        code = self.getInputFromPort('source')
        proc_code = urllib.unquote(str(code))
        cmd = 'def f(self,x):\n'
        for line in proc_code.split('\n'):
            cmd += '\t' + line + '\n'
        exec cmd
        
        if self.hasInputFromPort('xrange_min') and self.hasInputFromPort('xrange_max'):
            xmin = self.getInputFromPort('xrange_min')
            xmax = self.getInputFromPort('xrange_max')
            selection = (data.x >= xmin) & (data.x <= xmax)
            fw.fit(self,f, pars, data.y[selection], data.x[selection])
        else:
            fw.fit(self,f, pars, data.y, data.x)
        
        data.props['label'] = ''
        for i in range(0,len(pars)):
            p = raw_pars[i][0]
            data.props['label'] = data.props['label'] + p + ' = ' + str(pars[i].get()) + ', '
            data.props[p] = pars[i].get()
        
        data.x = np.linspace(min(data.x), max(data.x), 1000)
        data.y = f(self,data.x)
