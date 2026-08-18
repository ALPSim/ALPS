/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Mario Koenz <mkoenz@ethz.ch>                       *
 * Copyright (C) 2012 - 2014 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_NUMERIC_INF_HEADER
#define ALPS_NGS_NUMERIC_INF_HEADER

#include <limits>

namespace alps {
    namespace ngs {
        namespace numeric {

			template<typename T> struct inf {};

			template<> struct inf<double> {
        		operator double() const {
        			return std::numeric_limits<double>::infinity();
        		}
        	};

		}
	}
}

#endif