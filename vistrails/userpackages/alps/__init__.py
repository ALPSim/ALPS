# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Get ALPS at http://alps.comp-phys.org/
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

import dataset

identifier = 'org.comp-phys.alps'
version = '0.4.1'
name = 'ALPS'

configuration = ConfigurationObject(alpspath="/opt/alps/bin",toolpath="/opt/local/bin",mpirun="['mpirun','-np']",mpiprocs=0)

##############################################################################

def package_dependencies():
  return ['edu.utah.sci.vistrails.control_flow', 'edu.utah.sci.vistrails.matplotlib', 'edu.utah.sci.vistrails.spreadsheet']

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


