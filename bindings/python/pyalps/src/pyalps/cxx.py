# ****************************************************************************
#
# ALPS Project: Algorithms and Libraries for Physics Simulations
#
# ALPS Libraries
#
# Copyright (C) 2016 by Michele Dolfi <dolfim@phys.ethz.ch>
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
#
# ****************************************************************************


## Imports the compiled extensions from pyalps._ext and re-exports them, both
## as attributes of pyalps.cxx and as pyalps.cxx.<name> module entries.
##
## The historic second half of this file caught ImportError and retried the
## same names as top-level modules, for the old in-tree build where every .so
## sat side by side on PYTHONPATH. No installed layout can satisfy that spelling,
## so the fallback could only ever mask a genuine import failure.

# The compiled extensions, in the order pyalps.cxx exposes them.
_EXTENSIONS = (
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
)

import importlib
import sys

for _name in _EXTENSIONS:
    _module = importlib.import_module("." + _name, __package__ + "._ext")
    globals()[_name] = _module
    # Register under pyalps.cxx.<name> as well. This is not decoration: the
    # package's own modules (alea, alea_detail, hdf5, pytools, ngs) import
    # through `from .cxx.<name> import ...`, which is submodule syntax and
    # needs a sys.modules entry, as does any downstream code that spells it
    # the same way. Attribute access alone would not serve either.
    sys.modules["{}.cxx.{}".format(__package__, _name)] = _module

del _name, _module
