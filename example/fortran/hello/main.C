/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2011 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/parapack.h>
#include "alps/fortran/fortran_wrapper.h"

PARAPACK_SET_VERSION("my version");
PARAPACK_SET_COPYRIGHT("my copyright");
PARAPACK_REGISTER_WORKER(alps::fortran_wrapper, "hello fortran");
PARAPACK_REGISTER_EVALUATOR(alps::parapack::simple_evaluator, "hello fortran");

int main(int argc, char** argv) { return alps::parapack::start(argc, argv); }
