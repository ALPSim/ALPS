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

from .dataset import *
from .tools import *
from .pytools import *
from .floatwitherror import FloatWithError
from . import fit_wrapper
from . import cxx as cxx


# The extensions live in ``pyalps._ext`` in wheels, but Boost.Python-era
# installations also exposed the core modules directly below ``pyalps``.
# Register aliases instead of loading a second copy of an extension: nanobind
# has one process-wide type registry, and duplicate module instances would
# create subtly incompatible versions of the same C++ types.
for _extension_name in (
    "pyalea_c",
    "pymcdata_c",
    "pytools_c",
    "pyngsparams_c",
    "pyngshdf5_c",
    "pyngsbase_c",
    "pyngsobservable_c",
    "pyngsobservables_c",
    "pyngsresult_c",
    "pyngsresults_c",
    "pyngsapi_c",
    "pyngsrandom01_c",
    "pyngsaccumulator_c",
):
    _extension = getattr(cxx, _extension_name)
    globals()[_extension_name] = _extension
    sys.modules[__name__ + "." + _extension_name] = _extension

# Optional solver modules are present when the wheel was built from an ALPS
# checkout with application bindings enabled.
for _extension_name in ("maxent_c", "cthyb", "ctint"):
    try:
        _extension = __import__(
            __name__ + "._ext." + _extension_name, fromlist=[_extension_name]
        )
    except ImportError:
        continue
    globals()[_extension_name] = _extension
    sys.modules[__name__ + "." + _extension_name] = _extension

del _extension_name, _extension
