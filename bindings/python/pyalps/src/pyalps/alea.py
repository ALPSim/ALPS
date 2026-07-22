# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Olivier Parcollet
# 
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

""" 
This module contains classes for the evaluation of Monte Carlo measurements:
- RealObservable
- RealVectorObservable
- RealTimeSeriesObservable
- RealVectorTimeSeriesObservable

- MCScalarTimeseries
- MCVectorTimeseries
- MCScalarTimeseriesView
- MCVectorTimeseriesView
"""
from __future__ import print_function
from __future__ import absolute_import

from .cxx.pymcdata_c import *
from .cxx.pyalea_c import RealObservable, RealVectorObservable, RealTimeSeriesObservable, RealVectorTimeSeriesObservable
from .cxx.pyalea_c import MCScalarTimeseries, MCScalarTimeseriesView, MCVectorTimeseries, MCVectorTimeseries, ValueWithError, StdPairDouble, size, mean, variance, integrated_autocorrelation_time, running_mean, reverse_running_mean
from . import alea_detail as detail
import numpy
import pyalps.dataset


# TODO: __repr__ function for StdPairDouble
#def std_pair_double.__repr__():
#  return str(self.first) + ", " + str(self.second)

def make_dataset(MCTimeseries):
  "Makes a pyalps.dataset from a MCTimeseries object making it possible to plot one using pyalps.plot"
  d = pyalps.dataset.DataSet()
  d.y = MCTimeseries.timeseries()
  d.x = numpy.arange(1,len(d.y)+1,1)
  d.props['line'] = 'scatter'
  return d

def autocorrelation(timeseries, _distance = None, _limit = None):
  "Calculates the autocorrelation of a given timeseries. Can be invoked with _distance or _limit.\n\
timeseries: Any MCTimeseries or MCData object\n\
_distance: Calculates the autocorrelation until a specific length\n\
_limit: Calculates the autocorrelation until it has reached _limit of its initial value\n\
returns: MCTimeseries object with the autocorrelation"
  if _distance != None:
    return detail.autocorrelation_distance(timeseries, _distance)
  if _limit != None:
    return detail.autocorrelation_limit(timeseries, _limit)
  print("Usage: autocorrelation(timeseries, [_distance = XXX | _limit = XXX] )")

def cut_head(timeseries, _distance = None, _limit = None):
  "Creates a MCTimeseriesView object. Can be invoked with _distance or _limit.\n\
timeseries: Any MCTimeseries object\n\
_distance: Cuts the first _distance entries\n\
_limit: Cuts the front until the timeseries reaches _limit of its initial value\n\
returns: MCTimeseriesView object with the smaller timeseries\n\
  Note: does not copy the data, only creates a reference."
  if _distance != None:
    return detail.cut_head_distance(timeseries, _distance)
  if _limit != None:
    return detail.cut_head_limit(timeseries, _limit)
  print("Usage: cut_head(timeseries, [_distance = XXX | _limit = XXX] )")

def cut_tail(timeseries, _distance = None, _limit = None):
  "Creates a MCTimeseriesView object. Can be invoked with _distance or _limit.\n\
timeseries: Any MCTimeseries object\n\
_distance: Cuts the last _distance entries\n\
_limit: Cuts the tail until the timeseries only decays from its initial value to _limit of its initial value\n\
returns: MCTimeseriesView object with the smaller timeseries\n\
  Note: does not copy the data, only creates a reference."
  if _distance != None:
    return detail.cut_tail_distance(timeseries, _distance)
  if _limit != None:
    return detail.cut_tail_limit(timeseries, _limit)
  print("Usage: cut_head(timeseries, [_distance = XXX | _limit = XXX] )")

def exponential_autocorrelation_time(autocorrelation, _from = None, _to = None, _max = None, _min = None):
  "Fits a timeseries exponentially. Can be invoked with _from and _to or with _max and _min.\n\
autocorrelation: A MCTimeseries object\n\
_from & _to: fits the autocorrelation between _from and _to\n\
_max & _min: fits the autocorrelation between the values where it is at _max and where it is at _min from its initial value\n\
returns: StdPairDouble object FIT with the parameters of the fit\n\
  Note: The equation is FIT.first * exp(FIT.second * t)"
  if (_from != None and _to != None):
    return detail.exponential_autocorrelation_time_distance(autocorrelation, _from, _to)
  if (_max != None and _min != None):
    return detail.exponential_autocorrelation_time_limit(autocorrelation, _max, _min)
  print("Usage: exponential_autocorrelation_time(autocorrelation, [_from = XXX, _to = XXX | _max = XXX, _min = XXX] )")

binning = "binning"
uncorrelated = "uncorrelated"

def error(timeseries, selector = uncorrelated):
  "Calculates an error estimate of a timeseries. Can be invoked with 'uncorrelated' or 'binning'\n\
timeseries: Any MCTimeseries or MCData object\n\
uncorrelated (default): The error estimate when uncorrelated data is assumed\n\
binning: The error estimate is calculated using a binning analysis.\n\
returns: A float or numpyarray (depending on the dimension of the timeseries) with the error(s)"
  if selector == "binning":
    return detail.binning_error(timeseries)
  if selector == "uncorrelated":
    return detail.uncorrelated_error(timeseries)

