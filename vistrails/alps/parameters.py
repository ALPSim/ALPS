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

##############################################################################

class Parameters(Module):
   """ A dictionary of parameters """
   def update_parms(self,other):
       for k in other.parms.keys():
         if self.parms.has_key(k) and self.parms[k]!=other.parms[k]:
            raise ModuleError(self, "parameter " + k + " defined twice and the values differ")
         else: 
            self.parms[k] = other.parms[k]
       for k in other.defaults.keys():
         if self.defaults.has_key(k) and self.defaults[k]!=other.defaults[k]:
            raise ModuleError(self, "default for parameter " + k + " defined twice and the values differ")
         else: self.defaults[k] = other.defaults[k]
   def update_parms_unchecked(self,other):
       self.parms.update(other.parms)
       self.defaults.update(other.defaults)
   def update_parms_from_port(self,port_name):
       if self.hasInputFromPort(port_name):
         input_values = self.forceGetInputListFromPort(port_name)
         for q in input_values:
           self.update_parms(q)
   def set_checked(self,key_name,value):
       if type(value) == basic.File:
         self.set_checked(key_name,value.name)
       else:
         if self.parms.has_key(key_name) and self.parms[key_name] != value:
            raise ModuleError(self, "parameter " + key_name + " defined twice and the values differ")
         else:
           self.parms[key_name] = value
   def set_from_port(self,port_name):
       if self.hasInputFromPort(port_name):
          self.set_checked(port_name,self.getInputFromPort(port_name))
   def set_missing_from_defaults(self):
       for k in self.defaults.keys():
         if not self.parms.has_key(k):
           self.parms[k] = self.defaults[k]
   def readInputs(self):
       for port_name in self.inputPorts:
         if port_name != 'parms':
            self.set_from_port(port_name)
   def write(self,out):
     self.set_missing_from_defaults()
     out.write('{\n')
     for key in self.parms:
       out.write(key)
       out.write('=\"')
       out.write(self.parms[key])
       out.write('\"\n')
     out.write('}\n')
   def to_string(self):
     res = str(self.parms)
     if len(self.defaults)>0:
       res += '; defaults='
       res += str(self.defaults)
     return res
   def setOutput(self):
       self.update_parms_from_port('parms')
       self.setResult('value',self)
       self.setResult('value_as_string',self.to_string())
   def compute(self):
       self.parms = {}
       self.defaults = {}
       self.readInputs()
       self.setOutput()
   parms = {}
   defaults = {}

class FixedAndDefaultParameters(Parameters):
   def compute(self):
       self.parms = self.fixedparms.copy()
       self.defaults = self.fixeddefaults.copy()
       self.readInputs()
       self.setOutput()
   fixedparms={}
   fixeddefaults={}

class SystemParameters(Parameters):
   """ A dictionary of parameters defining a model and lattice and all their parameters 
   """

class UpdateParameters(Parameters):
   def compute(self):
       self.parms = {}
       self.defaults = {}
       if self.hasInputFromPort('parms'):
         input_values = self.forceGetInputListFromPort('parms')
         for p in input_values:
           self.update_parms_unchecked(p)
       if self.hasInputFromPort('updated_parms'):
         input_values = self.forceGetInputListFromPort('updated_parms')
         for p in input_values:
           self.update_parms_unchecked(p)
       self.setResult('value',self)
       self.setResult('value_as_string',str(self.parms))
   _input_ports = [('parms',[Parameters]),
                   ('updated_parms',[Parameters])]
                

class Parameter(Module):
   """ A module to define a single parameter """
   def compute(self):
           p = Parameters()
           p.parms = {self.getInputFromPort('name') : self.getInputFromPort('value')}
           self.setResult('value',p)
           self.setResult('value_as_string',str(p.parms))
   _input_ports = [('name',[basic.String]),
                   ('value',[basic.String])]
   _output_ports=[('value', [Parameters]),
                  ('value_as_string', [basic.String])]


class MonteCarloMeasurements(Parameters):
   """ a module to specify measurements for Monte Carlo simulations """
   def set_from_port(self,port_name):
       if self.hasInputFromPort(port_name):
          val  = self.getInputFromPort(port_name)
          if type(val) == basic.Boolean:
            self.set_checked('MEASURE['+port_name+']',val)
          else:
            self.set_checked(port_name,val)
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
       self.parms = {}
       if self.hasInputFromPort('T'):
          self.parms['T'] = self.getInputFromPort('T')
          if self.hasInputFromPort('beta'):
            raise ModuleError(self, "cannot define both T and beta")
       if self.hasInputFromPort('beta'):
          self.parms['beta'] = self.getInputFromPort('beta')
       self.setOutput()
   _input_ports = [('T',[(basic.String, 'the temperature')]),
                   ('beta',[(basic.String, 'the inverse temperature')])]



class ConservedQuantumnumbers(Parameters):
   """ defines conserved quantum numbers for exact diagonalization """
   _input_ports = [('CONSERVED_QUANTUMNUMBERS', [basic.String])]

class CustomMeasurements(Parameters): 
   """ definition of custom measurements for diagonalization and DMRG """

class CustomMeasurement(basic.Module):
   def compute(self):
           p = CustomMeasurements()
           key = 'MEASURE_'+ self.prefix+'['+str(self.getInputFromPort('name'))+']'
           p.parms = {key : self.getInputFromPort('definition')}
           self.setResult('value',p)
           self.setResult('value_as_string',str(p.parms))
   _input_ports = [('name',[(basic.String, 'the name of the measurement')]),
                   ('defintiion',[(basic.String, 'the operator to be measured')])]
   _output_ports=[('value', [CustomMeasurements]),
                  ('value_as_string', [basic.String])]

class AverageMeasurement(CustomMeasurement): 
   prefix = 'AVERAGE'

class LocalMeasurement(CustomMeasurement): 
   prefix = 'LOCAL'

class CorrelationsMeasurement(CustomMeasurement): 
   prefix = 'CORRELATIONS'
                   
class StructureFactorMeasurement(CustomMeasurement): 
   prefix = 'STRUCTURE_FACTOR'

class ParameterList(Module):
   """ a list of Parameters """
   def write(self,out):
     for p in self.list:
       p.write(out)
   def collectparmlist(self):
       if self.hasInputFromPort('parmlist'):
         input_values = self.forceGetInputListFromPort('parmlist')
         for p in input_values:
           self.list.extend(p.list)
   def collectparms(self):
       if self.hasInputFromPort('parms'):
         input_values = self.forceGetInputListFromPort('parms')
         for p in input_values:
           self.list.append(p)
   def to_string(self):
       strlist = []
       for p in self.list:
         strlist.append(p.parms)
       return str(strlist)
   def setOutput(self):
       self.setResult('value',self)
       self.setResult('value_as_string',self.to_string())   
   def compute(self):
       self.list=[]
       self.collectparmlist()
       self.collectparms()
       self.setOutput()
   list = []
   
class UpdateParameterList(ParameterList):
   def compute(self):
       self.list=[]
       self.collectparmlist()
       if self.hasInputFromPort('parms'):
         input_values = self.forceGetInputListFromPort('parms')
         for p in input_values:
           self.update_parms_unchecked(p)
       self.setOutput()

class IterateParameter(ParameterList):
   """ Iterate a parameter over a list of values """
   def compute(self):
       self.list=[]
       self.collectparmlist()
       self.collectparms()
       if self.hasInputFromPort('name') and self.hasInputFromPort('value_list'):
         key = self.getInputFromPort('name')
         value_list = self.getInputFromPort('value_list')
         input_values = self.forceGetInputListFromPort('parms')
         newlist = []
         for val in value_list:
           for p in self.list:
             newlist.append(SystemParameters())
             newlist[-1].parms=p.parms.copy()
             newlist[-1].parms[key]=str(copy.copy(val))
         self.list = newlist 
       self.setOutput()
   _input_ports = [('name',[basic.String]),
                   ('value_list',[ListOfElements])]

def register_parameters(type, ns="Parameters"):
  reg = core.modules.module_registry.get_module_registry()
  reg.add_module(type,namespace=ns)
  reg.add_output_port(type, "value", type)
  reg.add_output_port(type, "value_as_string", basic.String)

def register_abstract_parameters(type, ns="Parameters"):
  reg = core.modules.module_registry.get_module_registry()
  reg.add_module(type,namespace=ns,abstract=True)
  reg.add_output_port(type, "value", type)
  reg.add_output_port(type, "value_as_string", basic.String)
  
def initialize(): pass

def selfRegister():

  reg = core.modules.module_registry.get_module_registry()
    
  register_parameters(Parameters)
  reg.add_input_port(Parameters, "parms", Parameters)

  reg.add_module(Parameter,namespace="Parameters")
  
  reg.add_module(FixedAndDefaultParameters,namespace="Parameters",abstract=True)

  register_parameters(SystemParameters,"System")

  register_parameters(ParameterList)
  reg.add_input_port(ParameterList, "parmlist", ParameterList)
  reg.add_input_port(ParameterList, "parms", SystemParameters)
  
  register_parameters(UpdateParameters)
  register_parameters(UpdateParameterList)
  register_parameters(IterateParameter)
    
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
  

  
