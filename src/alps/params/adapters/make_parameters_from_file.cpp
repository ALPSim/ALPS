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
#include <alps/parameter.h>

#include "evaluated_parameters.hpp"

#include <boost/filesystem/fstream.hpp>

#include <sstream>
#include <stdexcept>
#include <string>

namespace alps {

    params make_parameters_from_file(boost::filesystem::path const & path) {
        boost::filesystem::ifstream ifs(path);
        std::ostringstream contents;
        contents << ifs.rdbuf();

        std::istringstream text(contents.str());
        try {
            return detail::make_parameters_evaluated(Parameters(text));
        } catch (std::runtime_error const & error) {
            // ParameterList holds one set per { } block, each with the keys set before it
            ParameterList sets;
            try {
                std::istringstream blocks(contents.str());
                sets.parse(blocks);
            } catch (std::runtime_error const &) {
                throw error;  // not blocks either
            }
            if (sets.size() != 1)
                throw std::runtime_error(path.string() + " has " + std::to_string(sets.size())
                    + " parameter sets in { } blocks, and a program reads one. "
                      "parameter2xml writes one XML job file per set.");
            return detail::make_parameters_evaluated(sets[0]);
        }
    }

}
