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
import platform

import dataset


##############################################################################

_subworkflows = [('MplXYPlotCell.xml', {'namespace': 'Tools'}),
                 ('ShowListOfPlots.xml', {'namespace': 'Dataset|Plot'}),
                 ('ShowMplPlot.xml', {'namespace': 'Dataset|Plot'}),
#                 ('ShowListOfPlots2.xml', {'namespace': 'DataSet|Plot', 'name':'ShowListOfPlots'}),
#                 ('ShowMplPlot2.xml', {'namespace': 'DataSet|Plot', 'name':'ShowMplPlot'}),
                 ('ShowListOfXMLFiles.xml', {'namespace': 'Tools'}),
                 ('ShowListOfHTMLFiles.xml', {'namespace': 'Tools'})]

def handle_module_upgrade_request(controller, module_id, pipeline):
   reg = get_module_registry()

   # format is {<old module name>: (<new_module_klass>, <remap_dictionary>}}
   # where remap_dictionary is {<remap_type>: <name_changes>}
   # and <name_changes> is a map from <old_name> to <new_name>
   module_remap = {'AlpsApplication': (applications.RunAlpsApplication,{}),
                   'AlpsEvaluate': (applications.AlpsEvaluate,{}),
                   'AppSpinMC': (applications.RunSpinMC,{}),
                   'AppLoop': (applications.RunLoop,{}),
                   'AppDirLoopSSE': (applications.RunDirLoopSSE,{}),
                   'AppWorm': (applications.RunWorm,{}),
                   'AppFullDiag': (applications.RunFullDiag,{}),
                   'AppSparseDiag': (applications.RunSparseDiag,{}),
                   'AppDMRG': (applications.RunDMRG,{}),
                   'AppTEBD': (applications.RunTEBD,{}),
                   'AppQWL': (applications.RunQWL,{}),
                   'EvaluateFullDiagT': (applications.EvaluateFullDiagVersusT,{}),
                   'EvaluateFullDiagH': (applications.EvaluateFullDiagVersusH,{}),
                   
                   'MonteCarloSimulation': (system.PrepareMonteCarlo,{}),
                   'DiagonalizationSimulation': (system.PrepareDiagonalization,{}),
                   'DMRGSimulation': (system.PrepareDMRG,{}),
                   'TEBDSimulation': (system.PrepareTEBD,{}),
                   'SimulationID': (system.SimulationName,{}),
                   'LatticeModel': (system.LatticeModel,{}),

                   'LatticeParameters': (lattices.Lattice,{}),
                   'square_lattice': (lattices.SquareLattice,{}),
                   'simple_cubic_lattice': (lattices.SimpleCubicLattice,{}),
                   'ladder': (lattices.LadderLattice,{}),
                   'open_ladder': (lattices.OpenLadderLattice,{}),
                   'chain_lattice': (lattices.ChainLattice,{}),
                   'open_chain_lattice': (lattices.OpenChainLattice,{}),
#                   'dimerized_chain_lattice': (lattices.DimerizedChainLattice,{}),
                   
                   'ModelParameters': (models.Model,{}),
                   'ClassicalSpinModel': (models.ClassicalSpinModel,{}),
                   'SpinModel': (models.SpinModel,{}),
                   'BosonHubbardModel': (models.BosonHubbardModel,{}),
                   'HardcoreBosonModel': (models.HardcoreBosonModel,{}),
                   
                   'CombineParameters': (parameters.ConcatenateParameters,{}),
                   'Parameter': (parameters.Parameter,{}),
                   'ConservedQuantumnumbers': (alpsparameters.ConservedQuantumNumbers,{}),
                   'SystemParameters': (alpsparameters.SystemParameters,{}),

                   'MakeParameterFile': (tools.WriteParameterFile,{}),
                   'MakeParameterXMLFiles': (tools.WriteInputFiles,{}),
                   'WriteInputFiles': (tools.WriteInputFiles,{}),
                   'GetRunFiles': (tools.GetCloneFiles,{}),
                   'GetResultFiles': (tools.GetResultFiles,{}),
                   'GetCloneFiles': (tools.GetCloneFiles,{}),
                   'XML2HTML': (tools.ConvertXML2HTML,{}),
                   'ConvertXML2HTML': (tools.ConvertXML2HTML,{}),
                   'GetSimulationInDir': (tools.GetJobFile,{}),
                   'OpenHTML': (alpscore.DisplayInBrowser,{}),
                   'TextFile': (alpscore.WriteTextFile,{}),

                   'GenerateDataSet': (dataset.PrepareDataSets,{}),
                   'LoadDataSet': (dataset.LoadDataSetsFromTextFile,{}),
                   'CustomLoader': (dataset.LoadCustomFile,{}),
                   'CollectXY': (dataset.CollectDataSets,{}),
                   'CollectDataSets': (dataset.CollectDataSets,{}),
                   'LoadProperties': (dataset.LoadAlpsProperties,{}),
                   'LoadAlpsHdf5': (dataset.LoadAlpsMeasurements,{}),
                   'LoadAlpsMeasurements': (dataset.LoadAlpsMeasurements,{}),
                   'LoadSpectrumHdf5': (dataset.LoadAlpsSpectra,{}),
                   'LoadBinningAnalysis': (dataset.LoadBinningAnalysis,{}),
                   'LoadAlpsDiagData': (dataset. LoadAlpsEigenstateMeasurements,{}),
                   'LoadAlpsEigenstateMeasurements': (dataset. LoadAlpsEigenstateMeasurements,{}),
                   'Transform': (dataset.TransformEachDataSet,{}),
                   'PlotDescriptor': (dataset.PreparePlot,{}),
                   'PreparePlot': (dataset.PreparePlot,{}),
                   'AxisDescriptor': (dataset.Axis,{}),
                   'Axis': (dataset.Axis,{}),
                   'LegendDescriptor': (dataset.Legend,{}),
                   'Legend': (dataset.Legend,{}),
                   'Convert2Text': (dataset.WriteTextFile,{}),
                   'Convert2Grace': (dataset.WriteGraceFile,{}),
                   'DisplayXMGRPlot': (plots.DisplayGracePlot,{}),
                   'GraceXYPlot': (dataset.WriteGraceFile,{}),
                   'MplXYPlot': (dataset.MplXYPlot,{}),
                   'Select': (dataset.Select,{}),
                   'And': (dataset.And,{}),
                   'Or': (dataset.Or,{}),
                   
                   'PolyFit': (dataset.DoPolynomialFit,{}),
                   'NonlinearFit': (dataset.DoNonlinearFit,{}),
                   'DoNonlinearFit': (dataset.DoNonlinearFit,{}),
                   
                   'SortByX': (dataset.SortEachDataSet,{}),
                   'SelectXRange': (dataset.RestrictXRange,{}),
                   'SetLabels': (dataset.SetLabels,{}),
                   'MakeScatter': (dataset.SetPlotStyle,{}),
                   'Selector': (dataset.Predicate,{}),
                   'PropertySelector': (dataset.PropertyPredicate,{}),
                   'PropertyRangeSelector': (dataset.PropertyRangePredicate,{}),
                   'ObservableSelector': (dataset.ObservablePredicate,{}),
                   'GroupBy': (dataset.GroupDataSets,{}),
                   'GroupDataSets': (dataset.GroupDataSets,{}),
                   'GroupedTransform': (dataset.TransformGroupedDataSets,{}),
                   'GenerateDataSet': (dataset.PrepareDataSets,{}),
                   'GenerateDataSet': (dataset.PrepareDataSets,{}),
                   
                   'CycleColors': (dataset.CycleColors,{}),
                   'CycleMarkers': (dataset.CycleMarkers,{}),
                   'Convert2XML': (tools.Convert2XML,{}),
                   'IterateValue': (parameters.IterateValue,{}),
                   'IterateParameter': (parameters.IterateParameter,{})
                   }


   new_remap = {}
   for name, (new_module, d) in module_remap.iteritems():
      new_remap[name] = [(None, '2.0.0', new_module, d)]

   # [(<start_version>, <end_version>, <new_module (None=same module, new version)>, <remap_dict>)]
   new_remap['ShowListOfHTMLFiles'] = [(None, '2.0.0', None, {})]
   new_remap['ShowListOfPlots'] = [(None, '2.0.0', None, {})]

   return UpgradeWorkflowHandler.remap_module(controller, module_id, pipeline,
                                             new_remap)


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
  
  alpscore.config = configuration
  
  dataset.initialize()


