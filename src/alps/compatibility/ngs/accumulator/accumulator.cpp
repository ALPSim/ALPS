/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2013 by Mario Koenz <mkoenz@ethz.ch>                       *
 *                              Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/accumulator.hpp>

namespace alps {
    namespace accumulator {

        #define EXTERNAL_FUNCTION(FUN)                          \
            result_wrapper FUN (result_wrapper const & arg) {   \
                return arg. FUN ();                             \
            }
            EXTERNAL_FUNCTION(sin)
            EXTERNAL_FUNCTION(cos)
            EXTERNAL_FUNCTION(tan)
            EXTERNAL_FUNCTION(sinh)
            EXTERNAL_FUNCTION(cosh)
            EXTERNAL_FUNCTION(tanh)
            EXTERNAL_FUNCTION(asin)
            EXTERNAL_FUNCTION(acos)
            EXTERNAL_FUNCTION(atan)
            EXTERNAL_FUNCTION(abs)
            EXTERNAL_FUNCTION(sqrt)
            EXTERNAL_FUNCTION(log)
            EXTERNAL_FUNCTION(sq)
            EXTERNAL_FUNCTION(cb)
            EXTERNAL_FUNCTION(cbrt)

        #undef EXTERNAL_FUNCTION
    }
}
