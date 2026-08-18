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

# The Boost.Python-era mpi_c extension is not part of the nanobind
# wheel build (no target builds it), so the old fallback chain
# (.cxx.mpi_c → mpi_c → boost.mpi) could never succeed anyway. Fail
# with an explanation instead of a misleading "No module named
# 'boost'".
raise ImportError(
    "pyalps.mpi is not available: the MPI bindings were not ported to the "
    "nanobind build of pyalps. Drive MPI-parallel simulations from C++, or "
    "use mpi4py for Python-side MPI communication."
)
