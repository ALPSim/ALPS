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

#include "ising.h"
#include <alps/parapack/parapack.h>
#include <alps/parapack/exchange.h>

int main(int argc, char** argv) { return alps::parapack::start(argc, argv); }

PARAPACK_SET_VERSION("ALPS/parapack example program: Classical Ising Model");
PARAPACK_REGISTER_ALGORITHM(ising_worker, "ising");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<ising_worker>, "ising; exchange");
PARAPACK_REGISTER_EVALUATOR(ising_evaluator, "ising");
