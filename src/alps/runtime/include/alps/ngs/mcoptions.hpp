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

#ifndef ALPS_NGS_MCOPTIONS_HPP
#define ALPS_NGS_MCOPTIONS_HPP

#include <alps/ngs/config.hpp>

#include <string>

namespace alps {

      class ALPS_DECL mcoptions {

        public:

            typedef enum { SINGLE, THREADED, MPI, HYBRID } execution_types;

            mcoptions(int argc, char* argv[]);

            bool valid;
            bool resume;
            std::size_t time_limit;
            std::string input_file;
            std::string output_file;
            std::string checkpoint_file;
            execution_types type;
    };
}

#endif
