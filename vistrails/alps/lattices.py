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

import core.modules.basic_modules
import core.modules.module_registry

import parameters

basic = core.modules.basic_modules

##############################################################################

class Lattice(parameters.FixedAndDefaultParameters):
   _input_ports = [('LATTICE',[basic.String]),
                   ('LATTICE_LIBRARY',[basic.File])]

class SquareLattice(Lattice):
   """ a square lattice """
   _input_ports = [('L',[(basic.String, 'the length')]),
                   ('W',[(basic.String, 'the width')])]
   fixed = {'LATTICE' : 'square lattice'}
   defaults = {'W':'L'}

class SimpleCubicLattice(Lattice):
   """ a simple cubic lattice """
   _input_ports = [('L',[(basic.String, 'the length')]),
                   ('W',[(basic.String, 'the width')]),
                   ('H',[(basic.String, 'the height')])]
   fixed = {'LATTICE' : 'simple cubic lattice'}
   defaults = {'W':'L', 'H':'L'}

class LadderLattice(Lattice):
   """ a LadderLattice lattice """
   _input_ports = [('L',[(basic.String, 'the length')]),
                   ('W',[(basic.String, 'the width')])]
   fixed = {'LATTICE' : 'ladder'}
   defaults = {'W':'2'}

class ChainLattice(Lattice):
   """ a chain lattice """
   _input_ports = [('L',[(basic.String, 'the length')])]
   fixed = {'LATTICE' : 'chain lattice'}

class OpenChainLattice(Lattice):
   """ an open chain lattice """
   _input_ports = [('L',[(basic.String, 'the length')])]
   fixed = {'LATTICE' : 'open chain lattice'}

class NNNChainLattice(Lattice):
    """ a chain lattice with nnn coupling """
    _input_ports = [('L',[(basic.String, 'the length')])]
    fixed = {'LATTICE' : 'nnn chain lattice'}

class DimerizedChainLattice(Lattice):
    """ a chain lattice with two different bonds """
    _input_ports = [('L',[(basic.String, 'the length')])]
    fixed = {'LATTICE' : 'dimerized chain lattice'}

class OpenLadderLattice(Lattice):
   """ an open LadderLattice lattice """
   _input_ports = [('L',[(basic.String, 'the length')]),
                   ('W',[(basic.String, 'the width')])]
   fixed = {'LATTICE' : 'open ladder'}
   defaults = {'W':'2'}


def register_lattice(type):
  reg = core.modules.module_registry.get_module_registry()
  reg.add_module(type,namespace="Lattices")
  reg.add_input_port(type,'LATTICE',[basic.String],True)

def initialize(): pass

def selfRegister():    

  reg = core.modules.module_registry.get_module_registry()

  reg.add_module(Lattice,namespace="Lattices")
  reg.add_output_port(Lattice, "value", Lattice)

  register_lattice(ChainLattice)
  register_lattice(OpenChainLattice)
  register_lattice(DimerizedChainLattice)
  register_lattice(LadderLattice)
  register_lattice(OpenLadderLattice)
  register_lattice(SquareLattice)
  register_lattice(SimpleCubicLattice)
  
 