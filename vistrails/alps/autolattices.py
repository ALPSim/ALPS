
# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Copyright (C) 2009 - 2010 by Bela Bauer <bauerb@itp.phys.ethz.ch>G
#
# Distributed under the Boost Software License, Version 1.0. (See accompany-
# ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
#
##############################################################################

import core.modules.basic_modules
import core.modules.module_registry

import parameters

basic = core.modules.basic_modules

##############################################################################

class Lattice(parameters.FixedAndDefaultParameters):
   """ a general lattice. Specify the lattice file in the LATTICE_LIBRARY input and the lattice name in LATTICE. LATTICE_LIBRARY defaults to the default ALPS lattices.xml file. """
   _input_ports = [('LATTICE',[basic.String]),
                   ('LATTICE_LIBRARY',[basic.File])]
class SquareLattice3x3(Lattice):
  """ automatically generated lattice: square lattice 3x3 """
  _input_ports = [
  ]
  fixed = {'lattice': 'square lattice 3x3'}

class SquareLattice4x4(Lattice):
  """ automatically generated lattice: square lattice 4x4 """
  _input_ports = [
  ]
  fixed = {'lattice': 'square lattice 4x4'}

class Dimer(Lattice):
  """ automatically generated lattice: dimer """
  _input_ports = [
  ]
  fixed = {'lattice': 'dimer'}

class Site(Lattice):
  """ automatically generated lattice: site """
  _input_ports = [
  ]
  fixed = {'lattice': 'site'}

class SimpleCubicLattice(Lattice):
  """ automatically generated lattice: simple cubic lattice """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')]),
    ('H',[(basic.String, '')])
  ]
  fixed = {'lattice': 'simple cubic lattice'}
  defaults = {'H': 'W', 'W': 'L'}

class SquareLattice(Lattice):
  """ automatically generated lattice: square lattice """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')])
  ]
  fixed = {'lattice': 'square lattice'}
  defaults = {'W': 'L'}

class OpenSquareLattice(Lattice):
  """ automatically generated lattice: open square lattice """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')])
  ]
  fixed = {'lattice': 'open square lattice'}
  defaults = {'W': 'L'}

class CoupledLadders(Lattice):
  """ automatically generated lattice: coupled ladders """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')])
  ]
  fixed = {'lattice': 'coupled ladders'}
  defaults = {'W': 'L'}

class TriangularLattice(Lattice):
  """ automatically generated lattice: triangular lattice """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')])
  ]
  fixed = {'lattice': 'triangular lattice'}
  defaults = {'W': 'L'}

class FrustratedSquareLattice(Lattice):
  """ automatically generated lattice: frustrated square lattice """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')])
  ]
  fixed = {'lattice': 'frustrated square lattice'}
  defaults = {'W': 'L'}

class ChainLattice(Lattice):
  """ automatically generated lattice: chain lattice """
  _input_ports = [
    ('L',[(basic.String, '')])
  ]
  fixed = {'lattice': 'chain lattice'}

class OpenChainLattice(Lattice):
  """ automatically generated lattice: open chain lattice """
  _input_ports = [
    ('L',[(basic.String, '')])
  ]
  fixed = {'lattice': 'open chain lattice'}

class NnnChainLattice(Lattice):
  """ automatically generated lattice: nnn chain lattice """
  _input_ports = [
    ('L',[(basic.String, '')])
  ]
  fixed = {'lattice': 'nnn chain lattice'}

class NnnOpenChainLattice(Lattice):
  """ automatically generated lattice: nnn open chain lattice """
  _input_ports = [
    ('L',[(basic.String, '')])
  ]
  fixed = {'lattice': 'nnn open chain lattice'}

class _2BandChainLattice(Lattice):
  """ automatically generated lattice: 2 band chain lattice """
  _input_ports = [
    ('L',[(basic.String, '')])
  ]
  fixed = {'lattice': '2 band chain lattice'}

class _2BandOpenChainLattice(Lattice):
  """ automatically generated lattice: 2 band open chain lattice """
  _input_ports = [
    ('L',[(basic.String, '')])
  ]
  fixed = {'lattice': '2 band open chain lattice'}

class AnisotropicSquareLattice(Lattice):
  """ automatically generated lattice: anisotropic square lattice """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')])
  ]
  fixed = {'lattice': 'anisotropic square lattice'}
  defaults = {'W': 'L'}

class Ladder(Lattice):
  """ automatically generated lattice: ladder """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')])
  ]
  fixed = {'lattice': 'ladder'}
  defaults = {'W': '2'}

class OpenLadder(Lattice):
  """ automatically generated lattice: open ladder """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')])
  ]
  fixed = {'lattice': 'open ladder'}
  defaults = {'W': '2'}

class InhomogeneousSquareLattice(Lattice):
  """ automatically generated lattice: inhomogeneous square lattice """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')])
  ]
  fixed = {'lattice': 'inhomogeneous square lattice'}
  defaults = {'W': 'L'}

class InhomogeneousSimpleCubicLattice(Lattice):
  """ automatically generated lattice: inhomogeneous simple cubic lattice """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')]),
    ('H',[(basic.String, '')])
  ]
  fixed = {'lattice': 'inhomogeneous simple cubic lattice'}
  defaults = {'H': 'W', 'W': 'L'}

class DepletedSquareLattice(Lattice):
  """ automatically generated lattice: depleted square lattice """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')])
  ]
  fixed = {'lattice': 'depleted square lattice'}
  defaults = {'W': 'L'}

class KagomeLattice(Lattice):
  """ automatically generated lattice: Kagome lattice """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')])
  ]
  fixed = {'lattice': 'Kagome lattice'}
  defaults = {'W': 'L'}

class HoneycombLattice(Lattice):
  """ automatically generated lattice: honeycomb lattice """
  _input_ports = [
    ('L',[(basic.String, '')]),
    ('W',[(basic.String, '')])
  ]
  fixed = {'lattice': 'honeycomb lattice'}
  defaults = {'W': 'L'}


def register_lattice(type):
  reg = core.modules.module_registry.get_module_registry()
  reg.add_module(type,namespace="Automatic lattices")
  reg.add_input_port(type,'LATTICE',[basic.String],True)

def initialize(): pass

def selfRegister():    

  reg = core.modules.module_registry.get_module_registry()
  
  reg.add_module(Lattice,namespace="Lattices")
  reg.add_output_port(Lattice, "value", Lattice)
  
  register_lattice(SquareLattice3x3)
  register_lattice(SquareLattice4x4)
  register_lattice(Dimer)
  register_lattice(Site)
  register_lattice(SimpleCubicLattice)
  register_lattice(SquareLattice)
  register_lattice(OpenSquareLattice)
  register_lattice(CoupledLadders)
  register_lattice(TriangularLattice)
  register_lattice(FrustratedSquareLattice)
  register_lattice(ChainLattice)
  register_lattice(OpenChainLattice)
  register_lattice(NnnChainLattice)
  register_lattice(NnnOpenChainLattice)
  register_lattice(_2BandChainLattice)
  register_lattice(_2BandOpenChainLattice)
  register_lattice(AnisotropicSquareLattice)
  register_lattice(Ladder)
  register_lattice(OpenLadder)
  register_lattice(InhomogeneousSquareLattice)
  register_lattice(InhomogeneousSimpleCubicLattice)
  register_lattice(DepletedSquareLattice)
  register_lattice(KagomeLattice)
  register_lattice(HoneycombLattice)
