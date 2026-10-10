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

#include "../single/ising.h"
#include <alps/parapack/exchange.h>
#ifdef ALPS_HAVE_MPI
#include "../multiple/ising.h"
#include <alps/parapack/exchange_multi.h>
#endif

PARAPACK_SET_VERSION("ALPS/parapack example program: exchange Monte Carlo");
PARAPACK_REGISTER_ALGORITHM(single_ising_worker, "ising");
PARAPACK_REGISTER_EVALUATOR(ising_evaluator, "ising");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<single_ising_worker>,
                         "ising; exchange");
#ifdef ALPS_HAVE_MPI
PARAPACK_REGISTER_PARALLEL_WORKER(alps::parapack::parallel_exchange_worker<single_ising_worker>,
                                  "ising; exchange");
PARAPACK_REGISTER_PARALLEL_WORKER(alps::parapack::multiple_parallel_exchange_worker<parallel_ising_worker>,
                                  "multiple parallel ising; exchange");
#endif
PARAPACK_REGISTER_EVALUATOR(ising_evaluator, "ising; exchange");
