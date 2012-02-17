#! /usr/bin/env python

from mcanalyze_tools import *

def calculate (obs):
  return alea.mean(obs)

impl_calculation("Mean", "mean/value", calculate)

