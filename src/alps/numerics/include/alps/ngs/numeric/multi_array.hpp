/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2014 by Jan Gukelberger <gukelberger@phys.ethz.ch>                *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_NUMERIC_MULTIARRAY_HEADER
#define ALPS_NGS_NUMERIC_MULTIARRAY_HEADER

#include <alps/multi_array/functions.hpp>

// Import multi_array functions into ngs::numeric namespace.
namespace alps {
    namespace ngs {
        namespace numeric {
            
            using alps::sin;
            using alps::cos;
            using alps::tan;
            using alps::sinh;
            using alps::cosh;
            using alps::tanh;
            using alps::asin;
            using alps::acos;
            using alps::atan;
            using alps::abs;
            using alps::sqrt;
            using alps::exp;
            using alps::log;
            using alps::fabs;

            using alps::sq;
            using alps::cb;
            using alps::cbrt;
            
            using alps::pow;
            using alps::sum;
        }
    }
}

#endif
