/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef PARAPACK_LOGGER_H
#define PARAPACK_LOGGER_H

#include <alps/parapack/process.h>
#include <alps/parapack/types.h>
#include <alps/utility/vmusage.hpp>
#include <iostream>
#include <string>

namespace alps {

struct logger {
  static std::string header();
  static std::string task(alps::tid_t tid);
  static std::string clone(alps::tid_t tid, alps::cid_t cid);
  static std::string group(alps::process_group g);
  static std::string group(alps::thread_group g);
  static std::string usage(alps::vmusage_type const& u);
};

} // namespace alps

#endif // PARAPACK_LOGGER_H
