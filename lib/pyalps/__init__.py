from __future__ import absolute_import
# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 1994-2009 by Bela Bauer <bauerb@phys.ethz.ch>
# 
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import sys
import os.path
if sys.platform == 'darwin' and not os.path.exists(os.path.expanduser('~/.matplotlib/matplotlibrc')):
    try:
        import matplotlib
        matplotlib.use('macosx')
    except ImportError:
        pass

from .dataset import *
from .tools import *
from .pytools import *
from .floatwitherror import FloatWithError
from . import fit_wrapper
