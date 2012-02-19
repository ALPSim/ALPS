# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2010 by Olivier Parcollet
# 
# This software is part of the ALPS libraries, published under the ALPS
# Library License; you can use, redistribute it and/or modify it under
# the terms of the license, either version 1 or (at your option) any later
# version.
#  
# You should have received a copy of the ALPS Library License along with
# the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
# available from http://alps.comp-phys.org/.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
# SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
# FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.
# 
# ****************************************************************************

""" 
This module contains classes for the evaluation of Monte Carlo measurements:
- RealObservable
- RealVectorObservable
- RealTimeSeriesObservable
- RealVectorTimeSeriesObservable
"""

from pyalea_c import *
from pymcdata_c import *

def autocorrelation(obs, _distance = None, _limit = None):
  if _distance != None:
    return autocorrelation_distance(obs, _distance)
  if _limit != None:
    return autocorrelation_limit(obs, _limit)
  print "Usage: autocorrelation(obs, [_distance = XXX | _limit = XXX] )"

def cut_head(obs, _distance = None, _limit = None):
  if _distance != None:
    return cut_head_distance(obs, _distance)
  if _limit != None:
    return cut_head_limit(obs, _limit)
  print "Usage: cut_head(obs, [_distance = XXX | _limit = XXX] )"

def cut_tail(obs, _distance = None, _limit = None):
  if _distance != None:
    return cut_tail_distance(obs, _distance)
  if _limit != None:
    return cut_tail_limit(obs, _limit)
  print "Usage: cut_head(obs, [_distance = XXX | _limit = XXX] )"

def exponential_autocorrelation_time(obs, _from = None, _to = None, _max = None, _min = None):
  if (_from != None and _to != None):
    return exponential_autocorrelation_time_distance(obs, _from, _to)
  if (_max != None and _min != None):
    return exponential_autocorrelation_time_limit(obs, _max, _min)
  print "Usage: exponential_autocorrelation_time(obs, [_from = XXX, _to = XXX | _max = XXX, _min = XXX] )"

binning = "binning"
uncorrelated = "uncorrelated"

def error(obs, selector = "uncorrelated"):
  if selector == "binning":
    return binning_error(obs)
  if selector == "uncorrelated":
    return uncorrelated_error(obs)

