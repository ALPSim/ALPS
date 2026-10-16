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

#ifndef ALPS_NGS_ALEA_ACCUMULATOR_ARGUMENTS_HPP
#define ALPS_NGS_ALEA_ACCUMULATOR_ARGUMENTS_HPP

#include <boost/parameter.hpp>

// = = = = N A M E D   P A R A M E T E R   D E F I N I T I O N = = = =

namespace alps {
    namespace accumulator  {

        BOOST_PARAMETER_NAME((bin_size, keywords) _bin_size)
        BOOST_PARAMETER_NAME((bin_num, keywords) _bin_num)
        BOOST_PARAMETER_NAME((weight_ref, keywords) _weight_ref)
        BOOST_PARAMETER_NAME((Weight, keywords) _Weight)

    } // end accumulator namespace
} // end alps namespace
#endif // ALPS_NGS_ALEA_ACCUMULATOR_ARGUMENTS_HPP
