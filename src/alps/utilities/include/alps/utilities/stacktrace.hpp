/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_STACKTRACE_HPP
#define ALPS_STACKTRACE_HPP

#include <alps/utilities_export.h>
#include <alps/ngs/config.hpp>
#include <alps/utilities/stringify.hpp>

#include <string>

// TODO: check for gcc and use __PRETTY_FUNCTION__

#define ALPS_STACKTRACE (                                                          \
      std::string("\nIn ") + __FILE__                                              \
    + " on " + ALPS_STRINGIFY(__LINE__)                                            \
    + " in " + __FUNCTION__ + "\n"                                          	   \
    + ::alps::stacktrace()                                                         \
)

namespace alps {

    ALPS_UTILITIES_DECL std::string stacktrace();

}

// Deprecated ALPS 2.x names, present only with the legacy layer (ALPS_LEGACY_API).
#if defined(__has_include)
#  if __has_include(<alps/legacy/ngs_stacktrace.hpp>)
#    include <alps/legacy/ngs_stacktrace.hpp>
#  endif
#endif

#endif
