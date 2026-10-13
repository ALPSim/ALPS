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

#include <alps/ngs/params_from_file.hpp>
#include <alps/parameter.h>

namespace alps {

    params params_from_file(boost::filesystem::path const & path) {
        params result;
        boost::filesystem::ifstream ifs(path);
        Parameters par(ifs);
        for (Parameters::const_iterator it = par.begin(); it != par.end(); ++it) {
            detail::paramvalue val(it->value());
            result[it->key()] = val;
        }
        return result;
    }

}
