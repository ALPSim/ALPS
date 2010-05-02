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

from core.configuration import ConfigurationObject
from core.upgradeworkflow import UpgradeWorkflowHandler
from core.modules.module_registry import get_module_registry

import alpscore
import parameters
import alpsparameters
import lattices
import models
import system
import applications
import plots
import tools
import evaluation
import platform

import dataset


##############################################################################

_subworkflows = [('MplXYPlotCell.xml', {'namespace': 'Tools'}),
                 ('ShowListOfPlots.xml', {'namespace': 'Dataset|Plot'}),
                 ('ShowMplPlot.xml', {'namespace': 'Dataset|Plot'}),
                 ('ShowListOfXMLFiles.xml', {'namespace': 'Tools'})]

def handle_module_upgrade_request(controller, module_id, pipeline):
   reg = get_module_registry()

   # format is {<old module name>: (<new_module_klass>, <remap_dictionary>}}
   # where remap_dictionary is {<remap_type>: <name_changes>}
   # and <name_changes> is a map from <old_name> to <new_name>
   module_remap = {'AlpsApplication': (applications.RunAlpsApplication,{}),
                   'AppSpinMC': (applications.RunSpinMC,{}),
                   'AppLoop': (applications.RunLoop,{}),
                   'AppDirLoopSSE': (applications.RunDirLoopSSE,{}),
                   'AppWorm': (applications.RunWorm,{}),
                   'AppFullDiag': (applications.RunFullDiag,{}),
                   'AppSparseDiag': (applications.RunSparseDiag,{}),
                   'AppDMRG': (applications.RunDMRG,{}),
                   'AppQWL': (applications.RunQWL,{}),
                   'EvaluateFullDiagT': (applications.EvaluateFullDiagVersusT,{}),
                   'EvaluateFullDiagH': (applications.EvaluateFullDiagVersusH,{}),
                   
                   'MonteCarloSimulation': (system.PrepareMonteCarlo,{}),
                   'DiagonalizationSimulation': (system.PrepareDiagonalization,{}),
                   'DMRGSimulation': (system.PrepareDMRG,{}),
                   'SimulationID': (system.SimulationName,{}),
                   'LatticeModel': (system.LatticeModel,{}),

                   'LatticeParameters': (lattices.Lattice,{}),
                   'square_lattice': (lattices.SquareLattice,{}),
                   'simple_cubic_lattice': (lattices.SimpleCubicLattice,{}),
                   'ladder': (lattices.LadderLattice,{}),
                   'open_ladder': (lattices.OpenLadderLattice,{}),
                   'chain_lattice': (lattices.ChainLattice,{}),
                   'open_chain_lattice': (lattices.OpenChainLattice,{}),
                   'dimerized_chain_lattice': (lattices.DimerizedChainLattice,{}),
                   
                   'ModelParameters': (models.Model,{}),
                   
                   'CombineParameters': (parameters.ConcatenateParameters,{}),
                   'ConservedQuantumnumbers': (alpsparameters.ConservedQuantumNumbers,{}),
                   'SystemParameters': (alpsparameters.SystemParameters,{}),

                   'MakeParameterFile': (tools.WriteParameterFile,{}),
                   'MakeParameterXMLFiles': (tools.WriteInputFiles,{}),
                   'GetRunFiles': (tools.GetCloneFiles,{}),
                   'XML2HTML': (tools.ConvertXML2HTML,{}),
                   'GetSimulationInDir': (tools.GetJobFile,{}),
                   'OpenHTML': (alpscore.DisplayInBrowser,{}),
                   'TextFile': (alpscore.WriteTextFile,{}),

                   'GenerateDataSet': (dataset.PrepareDataSets,{}),
                   'LoadDataSet': (dataset.LoadDataSetsFromTextFile,{}),
                   'CustomLoader': (dataset.LoadCustomFile,{}),
                   'CollectXY': (dataset.CollectDataSets,{}),
                   'LoadProperties': (dataset.LoadAlpsProperties,{}),
                   'LoadAlpsHdf5': (dataset.LoadAlpsMeasurements,{}),
                   'LoadSpectrumHdf5': (dataset.LoadAlpsSpectra,{}),
                   'LoadAlpsDiagData': (dataset. LoadAlpsEigenstateMeasurements,{}),
                   'Transform': (dataset.TransformEachDataSet,{}),
                   'PlotDescriptor': (dataset.PreparePlot,{}),
                   'AxisDescriptor': (dataset.Axis,{}),
                   'LegendDescriptor': (dataset.Legend,{}),
                   'Convert2Text': (dataset.WriteTextFile,{}),
                   'Convert2Grace': (dataset.WriteGraceFile,{}),
                   'DisplayXMGRPlot': (plots.DisplayGracePlot,{}),
                   'GraceXYPlot': (dataset.WriteGraceFile,{}),
                   
                   'PolyFit': (dataset.DoPolynomialFit,{}),
                   'NonlinearFit': (dataset.DoNonlinearFit,{}),
                   
                   'SortByX': (dataset.SortEachDataSet,{}),
                   'SelectXRange': (dataset.RestrictXRange,{}),
                   'SetLabels': (dataset.SetLabels,{}),
                   'MakeScatter': (dataset.SetPlotStyle,{}),
                   'Selector': (dataset.Predicate,{}),
                   'PropertySelector': (dataset.PropertyPredicate,{}),
                   'PropertyRangeSelector': (dataset.PropertyRangePredicate,{}),
                   'ObservableSelector': (dataset.ObservablePredicate,{}),
                   'GroupBy': (dataset.GroupDataSets,{}),
                   'GroupedTransform': (dataset.TransformGroupedDataSets,{}),
                   'GenerateDataSet': (dataset.PrepareDataSets,{}),
                   'GenerateDataSet': (dataset.PrepareDataSets,{}),
                   }


   old_module = pipeline.modules[module_id]
   print 'running module_upgrade_request', old_module.name
   if old_module.name in module_remap:
       print 'in_remap:', old_module.name
       remap = module_remap[old_module.name]
       new_descriptor = reg.get_descriptor(remap[0])
       try:
           function_remap = remap[1].get('function_remap', {})
           src_port_remap = remap[1].get('src_port_remap', {})
           dst_port_remap = remap[1].get('dst_port_remap', {})
           annotation_remap = remap[1].get('annotation_remap', {})
           action_list = \
               UpgradeWorkflowHandler.replace_module(controller, pipeline,
                                                     module_id, new_descriptor,
                                                     function_remap,
                                                     src_port_remap,
                                                     dst_port_remap,
                                                     annotation_remap)
       except Exception, e:
           import traceback
           traceback.print_exc()
           raise

       return action_list

   # otherwise, just try to automatic upgrade
   # attempt_automatic_upgrade
   return UpgradeWorkflowHandler.attempt_automatic_upgrade(controller, 
                                                           pipeline,
                                                           module_id)


def initialize():
  dataset.selfRegister()
  alpscore.selfRegister()  
  parameters.selfRegister()
  alpsparameters.selfRegister()
  lattices.selfRegister()
  models.selfRegister()
  applications.selfRegister()
  plots.selfRegister()
  tools.selfRegister()
  evaluation.selfRegister()
  
  alpscore.config = configuration
  
  dataset.initialize()


