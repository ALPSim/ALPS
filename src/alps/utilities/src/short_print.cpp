/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/short_print.hpp>

namespace alps { 
    namespace detail {
        std::ostream & operator<<(std::ostream & os, short_print_proxy<float> const & v) {
            std::streamsize precision = os.precision(v.precision);
            os << v.value;
            os.precision(precision);
            return os;
        }
    
        std::ostream & operator<<(std::ostream & os, short_print_proxy<double> const & v) {
            std::streamsize precision = os.precision(v.precision);
            os << v.value;
            os.precision(precision);
            return os;
        }

        std::ostream & operator<<(std::ostream & os, short_print_proxy<long double> const & v) {
            std::streamsize precision = os.precision(v.precision);
            os << v.value;
            os.precision(precision);
            return os;
        }
    } 
}
