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

#include "version.h"
#include <alps/utility/copyright.hpp>

std::string alps::parapack_copyright() {
  return "ALPS/parapack scheduler\n" \
    "  a Monte Carlo scheduler for multiple-level parallelization\n"    \
    "  copyright (c) 1997-" + alps::year() + " by Synge Todo <wistaria@comp-phys.org>\n";
}
