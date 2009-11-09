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
   fixed = {'LATTICE' : 'square lattice'}
   defaults = {'W':'L'}

class simple_cubic_lattice(LatticeParameters):
   """ a simple cubic lattice """
   _input_ports = [('L',[(basic.String, 'the length')]),
                   ('W',[(basic.String, 'the width')]),
                   ('H',[(basic.String, 'the height')])]
   fixed = {'LATTICE' : 'simple cubic lattice'}
   defaults = {'W':'L', 'H':'L'}

class ladder(LatticeParameters):
   """ a ladder lattice """
   _input_ports = [('L',[(basic.String, 'the length')]),
                   ('W',[(basic.String, 'the width')])]
   fixed = {'LATTICE' : 'ladder'}
   defaults = {'W':'2'}

class chain_lattice(LatticeParameters):
   """ a chain lattice """
   _input_ports = [('L',[(basic.String, 'the length')])]
   fixed = {'LATTICE' : 'chain lattice'}


def register_lattice(type):
  reg = core.modules.module_registry.get_module_registry()
  reg.add_module(type,namespace="Lattices")
  reg.add_input_port(type,'LATTICE',[basic.String],True)

def initialize(): pass

def selfRegister():    

  reg = core.modules.module_registry.get_module_registry()

  reg.add_module(LatticeParameters,namespace="Lattices")
  reg.add_output_port(LatticeParameters, "value", LatticeParameters)

  register_lattice(chain_lattice)
  register_lattice(ladder)
  register_lattice(square_lattice)
  register_lattice(simple_cubic_lattice)
  
 