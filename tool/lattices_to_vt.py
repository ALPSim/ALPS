#  Copyright Bela Bauer 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

import os, sys
from xml.etree import ElementTree

def replace_at(input, where, by):
    l = list(input)
    l[where] = by
    ret = ''
    for q in l:
        ret += q
    return ret

def camelify(s):
    if not s[0].isalpha():
        s = '_'+s
    s = s.strip()
    s = replace_at(s, 0, s[0].upper())
    while True:
        where = s.find(' ')
        if where == -1:
            break
        s = replace_at(s, where, '')
        s = replace_at(s, where, s[where].upper())
    return s

def write_lattice(name, inputs, fixed, defaults, outfile):
    outfile.write("class %s(Lattice):\n" % camelify(name))
    outfile.write('  """ automatically generated lattice: %s """\n' % name)
    outfile.write('  _input_ports = [\n')
    for i in range(len(inputs)):
        inp = inputs[i]
        outfile.write("    ('%s',[(basic.String, '')])" % inp)
        if i+1 == len(inputs):
            outfile.write("\n")
        else:
            outfile.write(",\n")
    outfile.write('  ]\n')
    if len(fixed) > 0:
        outfile.write('  fixed = %s\n' % fixed)
    if len(defaults) > 0:
        outfile.write('  defaults = %s\n' % defaults)
    outfile.write('\n')

def parse_latticegraphs(fn,outfile):
    root = ElementTree.parse(fn).getroot()
    names = []
    
    for lg in root.findall('LATTICEGRAPH'):
        name = lg.get('name')
        
        inputs = []
        fixed = {}
        defaults = {}
        
        fl = lg.find('FINITELATTICE')
        fixed['lattice'] = fl.find('LATTICE').get('ref')
        
        for extent in fl.findall('EXTENT'):
            size = extent.get('size')
            try:
                int(size)
                continue
            except ValueError:
                pass
            inputs.append(extent.get('size'))
        
        for param in fl.findall('PARAMETER'):
            defaults[param.get('name')] = param.get('default')
        
        # print name,inputs,fixed,defaults
        
        write_lattice(name, inputs, fixed, defaults, outfile)
        names.append(camelify(name))
    return names

if __name__ == '__main__':
    infilename = sys.argv[1]
    if sys.argv[2] == '-':
        outfile = sys.stdout
    else:
        outfile = open(sys.argv[2],'w')
    
    outfile.write('''
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
''')
    
    lattices = parse_latticegraphs(infilename, outfile)
    
    outfile.write('''
def register_lattice(type):
  reg = core.modules.module_registry.get_module_registry()
  reg.add_module(type,namespace="Automatic lattices")
  reg.add_input_port(type,'LATTICE',[basic.String],True)

def initialize(): pass

def selfRegister():    

  reg = core.modules.module_registry.get_module_registry()
  
  reg.add_module(Lattice,namespace="Lattices")
  reg.add_output_port(Lattice, "value", Lattice)
  
''')
    
    for lattice in lattices:
        outfile.write('  register_lattice(%s)\n' % lattice)
