/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef PARAPACK_UTIL_H
#define PARAPACK_UTIL_H

#include <alps/config.h>
#include <string>

namespace alps {

int hash(int n, int s = 826);

ALPS_DECL std::string id2string(int id, std::string const& pad = "_");

ALPS_DECL double parse_percentage(std::string const& str);

} // end namespace alps

#endif // PARAPACK_UTIL_H
