# VisTrails package for ALPS, Algorithms and Libraries for Physics Simulations
#
# Get ALPS at http://alps.comp-phys.org/
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

from parameters import Parameters

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
    _input_ports = [('REPRESENTATION',[(basic.String, 'the representation (SSE or PATHINTEGRAL)')])]

class ClassicalMonteCarloParameters(MonteCarloParameters): 
    """ A module to set the temperature """
    _input_ports = [('UPDATE',[(basic.String, 'the update method (local or cluster)')])]

class Temperature(Parameters):
    """ A module to set the temperature """
    def compute(self):
        if self.hasInputFromPort('T') and self.hasInputFromPort('beta'):
            raise ModuleError(self, "cannot define both T and beta")
        self.setOutput(self.readInputs(ParametersData({})))
    _input_ports = [('T',[(basic.String, 'the temperature')]),
                   ('beta',[(basic.String, 'the inverse temperature')])]


class ConservedQuantumnumbers(Parameters):
    """ defines conserved quantum numbers for exact diagonalization """
    _input_ports = [('CONSERVED_QUANTUMNUMBERS', [basic.String])]


class CustomMeasurements(Parameters): 
    """ definition of custom measurements for diagonalization and DMRG """


class CustomMeasurement(basic.Module):
    def compute(self):
        key = 'MEASURE_'+ self.prefix+'['+str(self.getInputFromPort('name'))+']'
        self.setResult('value',{key : self.getInputFromPort('definition')})
    _input_ports = [('name',[(basic.String, 'the name of the measurement')]),
                    ('defintiion',[(basic.String, 'the operator to be measured')])]
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

  register_parameters(SystemParameters,"System")

  register_parameters(MonteCarloMeasurements)
  register_parameters(MonteCarloParameters)
  reg.add_module(LoopMonteCarloParameters,namespace="Parameters")
  reg.add_module(ClassicalMonteCarloParameters,namespace="Parameters")
  register_parameters(Temperature)
  
  register_parameters(ConservedQuantumnumbers)
  register_abstract_parameters(CustomMeasurements)
  reg.add_module(CustomMeasurement,namespace="Parameters",abstract=True)
  reg.add_module(LocalMeasurement,namespace="Parameters")
  reg.add_module(AverageMeasurement,namespace="Parameters")
  reg.add_module(CorrelationsMeasurement,namespace="Parameters")
  reg.add_module(StructureFactorMeasurement,namespace="Parameters")
  

    
