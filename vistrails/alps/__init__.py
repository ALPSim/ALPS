# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Get ALPS at http://alps.comp-phys.org/
#
##############################################################################

from core.configuration import ConfigurationObject

import alpscore
import parameters
import lattices
import models
import system
import applications
import plots
import evaluation
import tools

identifier = 'org.comp-phys.alps'
version = '0.1.7'
name = 'ALPS'

configuration = ConfigurationObject(path=(None, str))

##############################################################################

def package_dependencies():
  return ['edu.utah.sci.vistrails.control_flow', 'edu.utah.sci.vistrails.matplotlib', 'edu.utah.sci.vistrails.spreadsheet']

def initialize():
  

  alpscore.selfRegister()  
  parameters.selfRegister()
  lattices.selfRegister()
  models.selfRegister()
  system.selfRegister()
  applications.selfRegister()
  plots.selfRegister()
  evaluation.selfRegister()
  tools.selfRegister()
  
  if configuration.check('path'):
    print "Have path", configuration.path
    alpscore.binpath = configuration.path

