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
from packages.alps.applications import RunAlpsApplication
import ising
import os

from packages.alps.parameters import Parameters 
basic = core.modules.basic_modules

identifier = 'org.comp-phys.alpsising'
version = '2.0.0'
name = 'Full ALPS Package tutorial'


def package_dependencies():
  return ['org.comp-phys.alps']

