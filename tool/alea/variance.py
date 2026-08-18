#! /usr/bin/env python
#/*****************************************************************************
#*
#* ALPS Project: Algorithms and Libraries for Physics Simulations
#*
#* Copyright (C) 2011-2012 by Lukas Gamper <gamperl@gmail.com>,
#*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
#*                            Maximilian Poprawe <poprawem@ethz.ch>
#*
#* ALPS Project: https://alps.comp-phys.org/
#* SPDX-License-Identifier: MIT
#*
#*****************************************************************************/

from mcanalyze_tools import *

def calculate (obs):
  return alea.variance(obs)

impl_calculation("Variance", "variance/value", calculate)

