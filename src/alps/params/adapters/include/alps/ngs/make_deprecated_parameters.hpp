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


#ifndef ALPS_NGS_MADE_DEPRECATED_PARAMETERS_HPP
#define ALPS_NGS_MADE_DEPRECATED_PARAMETERS_HPP

#include <alps/params.hpp>
#include <alps/config.h>
#include <alps/parameter.h>

namespace alps {

    /// The supplied parameters as an alps::Parameters, for code that uses it
    /** Each supplied parameter (params::supplied_text()) is passed on as
        its text, whether or not the program defined it. Defaults are not
        passed on: code reading the result applies its own. */
    ALPS_DECL Parameters make_deprecated_parameters(params const & arg);

}

#endif
