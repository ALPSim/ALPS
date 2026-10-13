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

#include <alps/ngs/make_deprecated_parameters.hpp>

#include <string>
#include <sstream>

namespace alps {

    Parameters make_deprecated_parameters(params const & arg) {
        Parameters par;
        for (params::const_iterator it = arg.begin(); it != arg.end(); ++it){
            std::stringstream s;
            s<<it->second;
            par.push_back(it->first,s.str());
        }
        return par;
    }
}
