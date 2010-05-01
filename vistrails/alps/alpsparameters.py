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
import copy
# import packages.controlflow
from packages.controlflow.list_module import ListOfElements
basic = core.modules.basic_modules

import parameters
from parameters import Parameters, ParametersData

##############################################################################

class SystemParameters(Parameters):
   """ A dictionary of parameters defining a model and lattice and all their parameters 
   """

class MonteCarloMeasurements(Parameters):
    """ a module to specify measurements for Monte Carlo simulations """
    def setFromPort(self,port_name,res):
        if self.hasInputFromPort(port_name):
          val  = self.getInputFromPort(port_name)
          print port_name," ",type(val)
          if type(val) == bool:
            res.set('MEASURE['+port_name+']',val)
          else:
            res.set(port_name,val)
        return res
            
    _input_ports = [('Correlations',[(basic.Boolean, 'the spin or density correlations')]),
                   ('Structure Factor',[(basic.Boolean, 'the spin or density structure factor')]),
                   ('Green Function',[(basic.Boolean, 'the Green function')]),
                   ('Site Type Density',[(basic.Boolean, 'the density at each type of site')]),
                   ('Bond Type Stiffness',[(basic.Boolean, 'the stiffness for each type of bond function')])
                   ]

class MonteCarloParameters(Parameters): 
    """ A module to set the Monte Carlo parameters """
    _input_ports = [('SWEEPS',[(basic.String, 'the number of sweeps for measurements')]),
                    ('THERMALIZATION',[(basic.String, 'the number of sweeps for thermalization')])
                    ]

class LoopMonteCarloParameters(MonteCarloParameters): 
    """ A module to set the temperature """
    def compute(self):
        res = self.readInputs(ParametersData({}))
        res.updateIfMissing(self.defaults)
        self.setResult('value',res.toBasic())
    _input_ports = [('ALGORITHM',[(basic.String, 'the algorithm to be used, default sse')])]
    defaults = { 'ALGORITHM':'sse' }

class QWLMonteCarloParameters(MonteCarloParameters): 
    """ A module to set the temperature """
    _input_ports = [('T_MIN',[basic.String]),
                    ('T_MAX',[basic.String]),
                    ('DELTA_T',[basic.String]),
                    ('CUTOFF',[basic.String]),
                    ('MEASURE_MAGNETIC_PROPERTIES',[basic.String]),
                    ('NUMBER_OF_WANG_LANDAU_STEPS',[basic.String],True),
                    ('USE_ZHOU_BHATT_METHOD',[basic.String],True),
                    ('FLATNESS_TRESHOLD',[basic.String],True),
                    ('BLOCK_SWEEPS',[basic.String],True),
                    ('INITIAL_MODIFICATION_FACTOR',[basic.String],True),
                    ('EXPANSION_ORDER_MINIMUM',[basic.String],True),
                    ('EXPANSION_ORDER_MAXIMUM',[basic.String],True),
                    ('START_STORING',[basic.String],True)
                    ]

class ClassicalMonteCarloParameters(MonteCarloParameters): 
    """ A module to set the temperature """
    _input_ports = [('UPDATE',[(basic.String, 'the update method (local or cluster)')])]


class DMRGParameters(Parameters): 
    """ A module to set the Monte Carlo parameters """
    _input_ports = [('SWEEPS',[(basic.String, 'the number of sweeps')]),
                    ('STATES',[(basic.String, 'the number of states in each sweep')]),
                    ('MAXSTATES',[(basic.String, 'the maximum number of states')]),
                    ('NUMBER_EIGENVALUES',[(basic.String, 'the number of eigenvalues to compute')])
                    ]


class Temperature(Parameters):
    """ A module to set the temperature """
    def compute(self):
        if self.hasInputFromPort('T') and self.hasInputFromPort('beta'):
            raise ModuleError(self, "cannot define both T and beta")
        self.setOutput(self.readInputs(ParametersData({})))
    _input_ports = [('T',[(basic.String, 'the temperature')]),
                   ('beta',[(basic.String, 'the inverse temperature')])]


class ConservedQuantumNumbers(Parameters):
    """ defines conserved quantum numbers for exact diagonalization """
    _input_ports = [('CONSERVED_QUANTUMNUMBERS', [basic.String])]


class CustomMeasurements(Parameters): 
    """ definition of custom measurements for diagonalization and DMRG """


class CustomMeasurement(basic.Module):
    def compute(self):
        print 'key:', self.getInputFromPort('name')
        print 'def:', self.getInputFromPort('definition')
        key = 'MEASURE_'+ self.prefix+'['+str(self.getInputFromPort('name'))+']'
        self.setResult('value',{key : self.getInputFromPort('definition')})
    _input_ports = [('name',[(basic.String, 'the name of the measurement')]),
                    ('definition',[(basic.String, 'the operator to be measured')])]
    _output_ports=[('value', [CustomMeasurements])]


class AverageMeasurement(CustomMeasurement): 
    prefix = 'AVERAGE'

class LocalMeasurement(CustomMeasurement): 
    prefix = 'LOCAL'

class CorrelationsMeasurement(CustomMeasurement): 
    prefix = 'CORRELATIONS'
                   
class StructureFactorMeasurement(CustomMeasurement): 
    prefix = 'STRUCTURE_FACTOR'



def register_parameters(type, ns="Parameters"):
    reg = core.modules.module_registry.get_module_registry()
    reg.add_module(type,namespace=ns)
    reg.add_output_port(type, "value", type)

def register_abstract_parameters(type, ns="Parameters"):
    reg = core.modules.module_registry.get_module_registry()
    reg.add_module(type,namespace=ns,abstract=True)
    reg.add_output_port(type, "value", type)
  
def initialize(): pass

def selfRegister():

  reg = core.modules.module_registry.get_module_registry()

  register_parameters(SystemParameters,"Parameters")

  register_parameters(MonteCarloMeasurements)
  register_parameters(MonteCarloParameters)
  register_parameters(DMRGParameters)
  reg.add_module(LoopMonteCarloParameters,namespace="Parameters")
  reg.add_module(ClassicalMonteCarloParameters,namespace="Parameters")
  reg.add_module(QWLMonteCarloParameters,namespace="Parameters")

  register_parameters(Temperature)
  
  register_parameters(ConservedQuantumNumbers)
  register_abstract_parameters(CustomMeasurements)
  reg.add_module(CustomMeasurement,namespace="Parameters",abstract=True)
  reg.add_module(LocalMeasurement,namespace="Parameters")
  reg.add_module(AverageMeasurement,namespace="Parameters")
  reg.add_module(CorrelationsMeasurement,namespace="Parameters")
  reg.add_module(StructureFactorMeasurement,namespace="Parameters")
  

    
