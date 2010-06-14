# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Copyright (C) 2009 - 2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
#                              Brigitte Surer <surerb@phys.ethz.ch>
#
# Distributed under the Boost Software License, Version 1.0. (See accompany-
# ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
#
##############################################################################

from core.modules.vistrails_module import Module
import core.modules.basic_modules
import os.path
import ising

basic = core.modules.basic_modules

identifier = 'org.comp-phys.ising'
version = '1.0.0'
name = 'ALPS Ising tutorial'


##############################################################################

def package_dependencies():
  return []


class IsingSimulation(Module):
    def compute(self): 
        result_file = self.interpreter.filePool.create_file(suffix='.h5')
        L = self.getInputFromPort('L')
        beta = self.getInputFromPort('beta')
        N = self.getInputFromPort('N')
        fname = self.getInputFromPort('fname')
        sim = ising.Simulation(beta,L)
        sim.run(N/2,N)
        sim.save(result_file.name)
        
        self.setResult('result_file', result_file)  

    _input_ports = [('L', [basic.Integer]),
                    ('beta', [basic.Float]), ('N', [basic.Integer]),('fname', [basic.String]) ]
    _output_ports = [('result_file', [basic.File])]

    
def initialize():
    reg = core.modules.module_registry.get_module_registry()
    reg.add_module(IsingSimulation,namespace="MyPackages")

