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

from alpscore import Dictionary
from core.modules.constant_configuration import StandardConstantWidget

##############################################################################

class CommonParametersFunctions:
    def set_checked(self,p,key,value):
        if isinstance(value,basic.File):
          value = value.name
        if p.has_key(key) and p[key] != value:
          raise ModuleError(self, "parameter " + key + " defined twice and the values differ")
        else:
          p[key] = value

    def update_checked(self,lhs,rhs):
        for k in rhs.keys():
          self.set_checked(lhs,k,rhs[k])
            

class ParametersData(CommonParametersFunctions):
    """ A dictionary of parameters """
    def __init__(self, d):
        self.parms = copy.deepcopy(d)
        
    def merge_one(self,p):
        self.update_checked(self.parms,p.parms)
      
    def copy(self):
        p = ParametersData(self.parms)
        return p
        
    def set(self,key_name,value):
        self.set_checked(self.parms,key_name,value)

    def updateIfMissing(self,defaults):
        for k in defaults.keys():
          if not self.parms.has_key(k):
            self.parms[k] = defaults[k]

    def write(self,out):
        out.write('{\n')
        for key in self.parms:
          out.write(key)
          out.write('=\"')
          out.write(str(self.parms[key]))
          out.write('\"\n')
        out.write('}\n')
        
    def update_unchecked(self,other):
        self.parms.update(other)

    def toBasic(self):
        return self.parms
        
    parms = {}

class ParameterListData(CommonParametersFunctions):
    """ a list of Parameters """
    def __init__(self, l):
        self.list = copy.deepcopy(l)

    def copy(self):
        p = ParameterListData(self.list)
        return p
        
    def merge_one(self,p):
        if isinstance(p,ParametersData):
          for q in self.list:
            self.update_checked(q,p.parms)
        else:
          newlist = []
          for p1 in self.list:
            for p2 in p.list:
              newlist.append(p1.copy())
              self.update_checked(newlist[-1],p2)
          self.list = newlist
        
    def set(self,key_name,value):
        for p in self.list:
          self.set_checked(p,key_name,value)
          
    def updateIfMissing(self,defaults):
        for k in defaults.keys():
          for p in self.list:
            if not p.has_key(k):
              p[k] = defaults[k]

    def write(self,out):
        for p in self.list:
          pd = ParametersData(p)
          pd.write(out)
       
    def toBasic(self):
        return self.list

    def update_unchecked(self,other):
        for p in self.list:
          p.update(other)
       
    list = []

def make_parameter_data(p):
    if isinstance(p,list):
      return ParameterListData(p)
    if isinstance(p,dict):
      return ParametersData(p)
    else:
      raise "unknown type passed to make_parameter_data"

class ParameterValue(basic.String):
    """ a class holding values for a parameter """
    @staticmethod
    def get_widget_class():
        return StandardConstantWidget

class ParameterValueList(list):
    """ a class holding values for a parameter """


class Parameters(Module):
    """ Simulation parameters """
    def merge(self,p1,p2):
        if isinstance(p1,ParametersData):
          res = p2
          res.merge_one(p1)
        else:
          res = p1
          res.merge_one(p2)
        return res
          
    def updateFromPort(self,port_name,res):
        if self.hasInputFromPort(port_name):
          input_values = self.forceGetInputListFromPort(port_name)
          for p in input_values:
            res = self.merge(res,make_parameter_data(p))
        return res
        
    def setFromPort(self,port_name,res):
        if self.hasInputFromPort(port_name):
          val = self.getInputFromPort(port_name)
          if (isinstance(val,ParameterValueList)):
            l = []
            for v in val:
              l.append({port_name:v})
            res = self.merge(res,ParameterListData(l))
          else:
            res.set(port_name,val)
          return res
    
    def readInputs(self,res):
        for port_name in self.inputPorts:
          if port_name == 'parms':
            res = self.updateFromPort('parms',res)
          else: 
            res = self.setFromPort(port_name,res)
        return res

#   def update_parms(self,other):
#       update_checked(self.parms,other.parms)
       
    def setOutput(self,res):
        self.setResult('value',res.toBasic())
                      
    def compute(self):
        self.setOutput(self.readInputs(ParametersData({})))

class MergeParameters(Module):
    """ merge several lists of parameters into one """
    def compute(self):
        l = []
        if self.hasInputFromPort('parms'):
          input_values = self.forceGetInputListFromPort('parms')
          for p in input_values:
            if isinstance(p,list):
              l.extend(p)
            else:
              l.append(p)
        if len(l)==1:
          self.setResult('value',l[0])
        else:
          self.setResult('value',l)
        
    _input_ports = [('parms',[Parameters])]
    _output_ports = [('value',[Parameters])]


class CombineParameters(Module,CommonParametersFunctions):
    """ join parameters from several lists, merging the parameters from each index """
    def compute(self):
        if self.hasInputFromPort('parms'):
          input_values = self.forceGetInputListFromPort('parms')
          l = input_values[0]
          for p in input_values[1:]:
            if len(p) != len(l):
              raise ModuleError(self, "Cannot join lists of different length")
            for i in range(len(l)):
              self.update_checked(l[i],p[i])
        if len(l)==1:
          self.setResult('value',l[0])
        else:
          self.setResult('value',l)
        
    _input_ports = [('parms',[Parameters])]
    _output_ports = [('value',[Parameters])]



class FixedAndDefaultParameters(Parameters):
    def setOutput(self,res):
        res.updateIfMissing(self.defaults)
        self.setResult('value',res.toBasic())
    def compute(self):
        res = self.readInputs(ParametersData(self.fixed))
        self.setOutput(res)
    fixed    = {}
    defaults = {}



class IterateValue(Module):
    """ Iterate a parameter over a list of values """
    def compute(self):
        self.setResult('value',ParameterValueList(self.getInputFromPort('value_list')))
    _input_ports = [('value_list',[ListOfElements])]
    _output_ports = [('value', [basic.String])]


class Parameter(Module):
    """ A module to define a single parameter """
    def compute(self):
        name = self.getInputFromPort('name')
        value = self.getInputFromPort('value')
        if (isinstance(value,ParameterValueList)):
          list = []
          for v in value:
            list.append({name:v})
          self.setResult('value',list)
        else:
          d = {name: value}
          self.setResult('value',d)
    _input_ports = [('name',[basic.String]),
                    ('value',[basic.String])]
    _output_ports=[('value', [Parameters])]


class IterateParameter(Module):
    """ Iterate a parameter over a list of values """
    def compute(self):
        if self.hasInputFromPort('name') and self.hasInputFromPort('value_list'):
          name = self.getInputFromPort('name')
          value_list = self.getInputFromPort('value_list')
          list = []
          for val in value_list:
            list.append({ name : val })
          self.list = list 
        self.setResult('value',list)
    _input_ports = [('name',[basic.String]),
                    ('value_list',[ListOfElements])]
    _output_ports = [('value', [Parameters])]


class UpdateParameters(Parameters):
   def compute(self):
       res=self.updateFromPort('parms',ParametersData({}))
       if self.hasInputFromPort('updated'):
         input_values = self.forceGetInputListFromPort('updated')
         for p in input_values:
           res.update_unchecked(p)
       self.setOutput(res)
   _input_ports = [('parms',[Parameters]),
                   ('updated',[Parameters])]


def register_parameters(type, ns="Parameters"):
  reg = core.modules.module_registry.get_module_registry()
  reg.add_module(type,namespace=ns)
  reg.add_output_port(type, "value", type)
  
def initialize(): pass

def selfRegister():

  reg = core.modules.module_registry.get_module_registry()

#  reg.add_module(ParametersData,namespace="Parameters")
#  reg.add_module(ParameterListData,namespace="Parameters")

  register_parameters(Parameters)
  reg.add_input_port(Parameters, "parms", Parameters)

  reg.add_module(ParameterValue,namespace="Parameters")
  reg.add_module(Parameter,namespace="Parameters")
  reg.add_module(FixedAndDefaultParameters,namespace="Parameters",abstract=True)
  reg.add_module(MergeParameters,namespace="Parameters")
  reg.add_module(CombineParameters,namespace="Parameters")
  reg.add_module(IterateParameter,namespace="Parameters")
  reg.add_module(IterateValue,namespace="Parameters")
  reg.add_module(UpdateParameters,namespace="Parameters")

