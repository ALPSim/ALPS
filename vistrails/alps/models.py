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

from core.modules.vistrails_module import Module, ModuleError, NotCacheable
import core.bundles
import core.modules.basic_modules
import core.modules.module_registry

import parameters

basic = core.modules.basic_modules

##############################################################################


class Model(parameters.FixedAndDefaultParameters): 
    """ A dictionary of parameters defining a model"""

class ClassicalSpinModel(Model):
    """ the classical spin models for ALPS spinmc """
    _input_ports = [('MODEL',[basic.String],True)]

class ClassicalIsingModel(ClassicalSpinModel):
    """ the classical Ising model for ALPS spinmc """
    fixed = {'MODEL':'Ising'}

class ClassicalXYModel(ClassicalSpinModel):
    """ the classical XY model for ALPS spinmc """
    fixed = {'MODEL':'XY'}

class ClassicalHeisenbergModel(ClassicalSpinModel):
    """ the classical Heisenberg model for ALPS spinmc """
    fixed = {'MODEL':'Heisenberg'}

                    
class SpinModel(Model):
   fixed = {'MODEL'   : 'spin'}
   defaults =      {'J0'    : '0',
                    'J'     : 'J0',
                    'Jxy'   : 'J',
                    'Jz'    : 'J',
                    'Jxy0'  : 'Jxy',
                    'Jz0'   : 'Jz',
                    'J1'    : '0',
                    "J'"    : 'J1',
                    "Jz'"   : "J'",
                    "Jxy'"  : "J'",
                    'Jz1'   : "Jz'",
                    'Jxy1'  : "Jxy'",
                    'h'     : '0',
                    'Gamma' : '0',
                    'D'     : '0',
                    'h#'    : 'h',
                    'Gamma#': 'Gamma',
                    'D#'    : 'D',
                    'J#'    : '0',
                    'Jz#'   : 'J#',
                    'Jxy#'  : 'J#'
                 }

class BosonHubbardModel(Model):
   fixed = {'MODEL'   : 'boson Hubbard'}
   defaults =      {'mu'    : '0',
                    't'     : '1',
                    'V'     : '0',
                    "t'"    : '0',
                    "V''"   : '0',
                    'U'     : '0',
                    't0'    : 't',
                    't1'    : "t'",
                    'V0'    : '0',
                    'V1'    : "V'",
                    'mu#'   : 'mu',
                    'U#'    : 'U',
                    't#'    : '0',
                    'V#'    : '0'
                 }

class HardcoreBosonModel(Model):
   fixed = {'MODEL'   : 'hardcore boson'}
   defaults =      {'mu'    : '0',
                    't'     : '1',
                    'V'     : '0',
                    "t'"    : '0',
                    "V''"   : '0',
                    't0'    : 't',
                    't1'    : "t'",
                    'V0'    : '0',
                    'V1'    : "V'",
                    'mu#'   : 'mu',
                    't#'    : '0',
                    'V#'    : '0'
                 }

def register_model(type):
   reg = core.modules.module_registry.get_module_registry()
   reg.add_module(type,namespace="Models")
   reg.add_input_port(type,'MODEL_LIBRARY',[basic.File])

def register_parameters(type, ns="Models"):
   reg = core.modules.module_registry.get_module_registry()
   reg.add_module(type,namespace=ns)
   reg.add_output_port(type, "value", type)
  
def initialize(): pass

def selfRegister():

   reg = core.modules.module_registry.get_module_registry()
  
   register_parameters(Model)
   
   reg.add_module(ClassicalSpinModel,namespace="Models",abstract=True)
   reg.add_module(ClassicalIsingModel,namespace="Models")
   reg.add_module(ClassicalXYModel,namespace="Models")
   reg.add_module(ClassicalHeisenbergModel,namespace="Models")
   
   register_model(SpinModel)
   register_model(BosonHubbardModel)
   register_model(HardcoreBosonModel)
  
