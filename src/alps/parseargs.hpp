/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_SCHEDULER_PARSEARGS_HPP
#define ALPS_NGS_SCHEDULER_PARSEARGS_HPP

#include <alps/ngs/config.hpp>

#include <string>

namespace alps {

	struct ALPS_DECL parseargs {
	    parseargs(int argc, char *argv[]);

	    bool resume;

	    std::size_t timelimit;
	    std::size_t tmin;
	    std::size_t tmax;

	    std::string input_file;
	    std::string output_file;
	};

}

#endif
