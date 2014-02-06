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

from vistrails.core.modules.vistrails_module import Module
import vistrails.core.modules.basic_modules
import ising

basic = vistrails.core.modules.basic_modules

class IsingSimulation(Module):
    def compute(self): 
        result_file = self.interpreter.filePool.create_file(suffix='.h5')
        L = self.getInputFromPort('L')
        beta = self.getInputFromPort('beta')
        N = self.getInputFromPort('N')
        sim = ising.Simulation(beta,L)
        sim.run(N/2,N)
        sim.save(result_file.name)
        
        self.setResult('result_file', result_file)  

    _input_ports = [('L', [basic.Integer]),
                    ('beta', [basic.Float]), ('N', [basic.Integer]) ]
    _output_ports = [('result_file', [basic.File])]

    
_modules = [IsingSimulation]
