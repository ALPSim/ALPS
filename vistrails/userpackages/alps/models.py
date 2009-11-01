# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Get ALPS at http://alps.comp-phys.org/
#
##############################################################################

from core.modules.vistrails_module import Module, ModuleError, NotCacheable
import core.bundles
import core.modules.basic_modules
import core.modules.module_registry

import parameters

basic = core.modules.basic_modules

##############################################################################


class ModelParameters(parameters.FixedAndDefaultParameters): 
    """ A dictionary of parameters defining a model"""
    _input_ports = [('MODEL',[basic.String]),
                    ('MODEL_LIBRARY',[basic.File])]

class ClassicalSpinModel(ModelParameters):
    """ the classical spin models for ALPS spinmc """

class SpinModel(ModelParameters):
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

class BosonHubbardModel(ModelParameters):
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

def register_model(type):
   reg = core.modules.module_registry.get_module_registry()
   reg.add_module(type,namespace="Models")
   reg.add_input_port(type,'MODEL',[basic.String],True)

def register_parameters(type, ns="Models"):
   reg = core.modules.module_registry.get_module_registry()
   reg.add_module(type,namespace=ns)
   reg.add_output_port(type, "value", type)
  
def initialize(): pass

def selfRegister():

   reg = core.modules.module_registry.get_module_registry()
  
   register_parameters(ModelParameters)
   reg.add_module(ClassicalSpinModel,namespace="Models")
   register_model(SpinModel)
   register_model(BosonHubbardModel)
  
