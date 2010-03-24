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

import alpscore
import parameters
import alpsparameters
import lattices
import models
import system
import applications
import plots
import tools
import evaluation
import platform

import dataset


##############################################################################

_subworkflows = [('MplXYPlotCell.xml', {'namespace': 'Tools'}) ]
#_subworkflows = ['MplXYPlotCell.xml']


def initialize():
  dataset.selfRegister()
  alpscore.selfRegister()  
  parameters.selfRegister()
  alpsparameters.selfRegister()
  lattices.selfRegister()
  models.selfRegister()
  system.selfRegister()
  applications.selfRegister()
  plots.selfRegister()
  tools.selfRegister()
  evaluation.selfRegister()
  
  alpscore.config = configuration
  
  dataset.initialize()


