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

#include <alps/ngs/make_parameters_from_file.hpp>
#include <alps/ngs/make_parameters_from_deprecated.hpp>
#include <alps/parameter.h>

#include <boost/filesystem/fstream.hpp>

namespace alps {

    params make_parameters_from_file(boost::filesystem::path const & path) {
        boost::filesystem::ifstream ifs(path);
        return make_parameters_from_deprecated(Parameters(ifs));
    }

}
