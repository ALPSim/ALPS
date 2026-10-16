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

#ifndef ALPS_STRINGIFY_HPP
#define ALPS_STRINGIFY_HPP

#define ALPS_STRINGIFY(arg) ALPS_STRINGIFY_HELPER(arg)

#define ALPS_STRINGIFY_HELPER(arg) #arg

// Deprecated ALPS 2.x names, present only with the legacy layer (ALPS_LEGACY_API).
#if defined(__has_include)
#  if __has_include(<alps/legacy/ngs_stringify.hpp>)
#    include <alps/legacy/ngs_stringify.hpp>
#  endif
#endif

#endif
