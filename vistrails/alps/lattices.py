# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Get ALPS at http://alps.comp-phys.org/
#
##############################################################################

import core.modules.basic_modules
import core.modules.module_registry

import parameters

basic = core.modules.basic_modules

##############################################################################

class LatticeParameters(parameters.FixedAndDefaultParameters):
   _input_ports = [('LATTICE',[basic.String]),
                   ('LATTICE_LIBRARY',[basic.File])]

class square_lattice(LatticeParameters):
   """ a square lattice """
   _input_ports = [('L',[(basic.String, 'the length')]),
                   ('W',[(basic.String, 'the width')])]
   fixedparms = {'LATTICE' : 'square lattice'}
   fixeddefaults = {'W':'L'}

class ladder(LatticeParameters):
   """ a ladder lattice """
   _input_ports = [('L',[(basic.String, 'the length')]),
                   ('W',[(basic.String, 'the width')])]
   fixedparms = {'LATTICE' : 'ladder'}
   fixeddefaults = {'W':'2'}

class chain_lattice(LatticeParameters):
   """ a chain lattice """
   _input_ports = [('L',[(basic.String, 'the length')])]
   fixedparms = {'LATTICE' : 'chain lattice'}

  
def register_lattice(type):
  reg = core.modules.module_registry.get_module_registry()
  reg.add_module(type,namespace="Lattices")
  reg.add_input_port(type,'LATTICE',[basic.String],True)

def initialize(): pass

def selfRegister():    

  reg = core.modules.module_registry.get_module_registry()

  reg.add_module(LatticeParameters,namespace="Lattices")
  reg.add_output_port(LatticeParameters, "value", LatticeParameters)
  reg.add_output_port(LatticeParameters, "value_as_string", basic.String)

  register_lattice(square_lattice)
  register_lattice(chain_lattice)
  register_lattice(ladder)
  
 