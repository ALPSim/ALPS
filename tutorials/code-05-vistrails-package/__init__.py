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

import core.modules.basic_modules
from packages.alps.applications import RunAlpsApplication

identifier = 'org.comp-phys.spinmc'
version = '1.0.0'
name = 'ALPS Package tutorial'


##############################################################################

def package_dependencies():
  return ['org.comp-phys.alps']


class RunSpinMonteCarlo(RunAlpsApplication):
    """Runs the spinmc classical Monte Carlo application """
    appname = 'spinmc'

    
def initialize():
    reg = core.modules.module_registry.get_module_registry()
    reg.add_module(RunSpinMonteCarlo,namespace="MyPackages")

