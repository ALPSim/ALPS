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
#include <alps/parapack/temperature_scan.h>

PARAPACK_SET_VERSION("ALPS/parapack example program: Metropolis algorithm");
PARAPACK_REGISTER_ALGORITHM(single_ising_worker, "ising");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::temperature_scan_adaptor<single_ising_worker>,
                         "ising; temperature scan");
PARAPACK_REGISTER_EVALUATOR(ising_evaluator, "ising");
