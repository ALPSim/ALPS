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
from packages.alps.applications import RunAlpsApplication
import ising
import os

from packages.alps.parameters import Parameters 
basic = vistrails.core.modules.basic_modules

class IsingAlpsSimulation(Module):
    def compute(self): 
        result_file = self.interpreter.filePool.create_file().name
        os.unlink(result_file)
        os.mkdir(result_file)
        dir=basic.Directory()
        dir.name=result_file
        list_of_parms = self.getInputFromPort('parm')
        for entry in list_of_parms:
            L = int(entry['L'])
            beta = float(entry['BETA'])
            N = int(entry['N'])
            sim = ising.Simulation(beta,L)
            sim.run(N/2,N)
            sim.save(os.path.join(result_file,"foobar.L%s_beta%.4s.h5"%(L,beta)))
        
        self.setResult('dir',dir)  

    _input_ports = [('parm', [Parameters])]
    _output_ports = [('dir', [basic.Directory])]

    
_modules = [IsingAlpsSimulation]