/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010-2012 by Haruhiko Matsuo <halm@looper.t.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Return virtual memory usage information (works only for Linux)
  
#ifndef ALPS_UTILITY_VMUSAGE_HPP
#define ALPS_UTILITY_VMUSAGE_HPP

#include <alps/config.h>
#include <map>
#include <string>


namespace alps {

typedef std::map<std::string, unsigned long> vmusage_type;

ALPS_DECL vmusage_type vmusage(int pid = -1);

} // end namespace alps

#endif // ALPS_UTILITY_VMUSAGE_HPP
