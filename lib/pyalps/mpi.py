# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 2012 by Matthias Troyer
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import sys
if sys.platform == 'linux2':
    import DLFCN as dl
    flags = sys.getdlopenflags()
    sys.setdlopenflags(dl.RTLD_NOW|dl.RTLD_GLOBAL)
    try:
        try:
            from .cxx.mpi_c import *
        except ImportError:
            from mpi_c import *
    except ImportError:
        from boost.mpi import *
    sys.setdlopenflags(flags)
else:
    try:
        try:
            from .cxx.mpi_c import *
        except ImportError:
            from mpi_c import *
    except ImportError:
        from boost.mpi import *