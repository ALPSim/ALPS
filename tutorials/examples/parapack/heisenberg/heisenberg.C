/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2012 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "heisenberg.h"
#include <alps/parapack/parapack.h>
#include <alps/parapack/exchange.h>

int main(int argc, char** argv) { return alps::parapack::start(argc, argv); }

PARAPACK_SET_VERSION("ALPS/parapack example program: Classical Heisenberg Model");
PARAPACK_REGISTER_ALGORITHM(heisenberg_worker, "heisenberg");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<heisenberg_worker>, "heisenberg; exchange");
PARAPACK_REGISTER_EVALUATOR(heisenberg_evaluator, "heisenberg");
