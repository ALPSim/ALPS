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

#ifndef ALPS_NGS_STACKTRACE_HPP
#define ALPS_NGS_STACKTRACE_HPP

#include <alps/ngs/config.hpp>
#include <alps/ngs/stringify.hpp>

#include <string>

// TODO: check for gcc and use __PRETTY_FUNCTION__

#define ALPS_STACKTRACE (                                                          \
      std::string("\nIn ") + __FILE__                                              \
    + " on " + ALPS_NGS_STRINGIFY(__LINE__)                                        \
    + " in " + __FUNCTION__ + "\n"                                          	   \
    + ::alps::ngs::stacktrace()                                                    \
)

namespace alps {
    namespace ngs {

        ALPS_DECL std::string stacktrace();

    }
}

#endif
