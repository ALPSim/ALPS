 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 #                                                                                 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2013 by Matthias Troyer <gamperl@gmail.com>                #
 #                                                                                 #
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
 #                                                                                 #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

try:
    from .maxent_c import *
except ImportError:
    from maxent_c import *
    