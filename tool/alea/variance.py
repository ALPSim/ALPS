#! /usr/bin/env python

from mcanalyze_tools import *

def calculate (obs):
  return alea.variance(obs)

impl_calculation("Variance", "variance/value", calculate)

