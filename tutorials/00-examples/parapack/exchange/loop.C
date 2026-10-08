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

#include "../loop/loop.h"
#include <alps/parapack/exchange.h>

PARAPACK_SET_VERSION("ALPS/parapack example program: exchange Monte Carlo");
PARAPACK_REGISTER_ALGORITHM(qmc_worker, "loop");
PARAPACK_REGISTER_EVALUATOR(alps::parapack::simple_evaluator, "loop");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<qmc_worker>, "loop; exchange");
#ifdef ALPS_HAVE_MPI
PARAPACK_REGISTER_PARALLEL_WORKER(alps::parapack::parallel_exchange_worker<qmc_worker>,
                                  "loop exchange");
#endif
PARAPACK_REGISTER_EVALUATOR(alps::parapack::simple_evaluator, "loop; exchange");
