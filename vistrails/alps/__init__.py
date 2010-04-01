# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Copyright (C) 2009 - 2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
#                              Synge Todo <wistaria@comp-phys.org>
#
# Distributed under the Boost Software License, Version 1.0. (See accompany-
# ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
#
##############################################################################

from core.configuration import ConfigurationObject
import platform

identifier = 'org.comp-phys.alps'
version = '0.5.1'
name = 'ALPS'


##############################################################################

def package_dependencies():
  return ['edu.utah.sci.vistrails.control_flow', 'edu.utah.sci.vistrails.matplotlib', 'edu.utah.sci.vistrails.spreadsheet', 'edu.utah.sci.vistrails.vtlcreator']


if platform.system()=='Windows':
  configuration = ConfigurationObject(alpspath="C:\\Program Files\\ALPS\\bin",toolpath="C:\\Program Files\\ALPS\\bin",mpirun="",mpiprocs=0)
else:
  configuration = ConfigurationObject(alpspath="/opt/alps/bin",toolpath="/opt/local/bin",mpirun="['mpirun','-np']",mpiprocs=0)

