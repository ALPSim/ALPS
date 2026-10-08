/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "loop.h"

PARAPACK_SET_VERSION("ALPS/parapack example program: loop cluster algorithm");
PARAPACK_REGISTER_ALGORITHM(qmc_worker, "loop");
PARAPACK_REGISTER_EVALUATOR(alps::parapack::simple_evaluator, "loop");
