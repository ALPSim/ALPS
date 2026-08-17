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

#include <alps/ngs/detail/paramproxy.hpp>

#include <alps/numeric/vector_functions.hpp>

#include <stdexcept>

namespace alps {
    namespace detail {

        void paramproxy::save(hdf5::archive & ar) const {
            if (!defined)
                throw std::runtime_error("No reference to parameter '" + key + "' available" + ALPS_STACKTRACE);
            ar[""] << (!!value ? *value : getter());
        }

        void paramproxy::load(hdf5::archive & ar) {
            if (!defined || !!value)
                throw std::runtime_error("No reference to parameter '" + key + "' available" + ALPS_STACKTRACE);
            if (!!value) {
                detail::paramvalue value;
                ar[""] >> value;
                setter(value);
            } else
                ar[""] >> *value;
        }

        void paramproxy::print(std::ostream & os) const {
			if (!defined)
				throw std::runtime_error("No parameter '" + key + "' available" + ALPS_STACKTRACE);
			os << (!value ? getter() : *value);
        }

        std::ostream & operator<<(std::ostream & os, paramproxy const & v) {
			v.print(os);
            return os;
		}

        #define ALPS_NGS_PARAMPROXY_ADD_OPERATOR_IMPL(T)                                \
            T operator+(paramproxy const & p, T s) {                                    \
                using boost::numeric::operators::operator+=;                            \
                return s += p.cast< T >();                                              \
            }                                                                           \
            T operator+(T s, paramproxy const & p) {                                    \
                using boost::numeric::operators::operator+=;                            \
                return s += p.cast< T >();                                              \
            }
        ALPS_NGS_FOREACH_PARAMETERVALUE_TYPE(ALPS_NGS_PARAMPROXY_ADD_OPERATOR_IMPL)
        #undef ALPS_NGS_PARAMPROXY_ADD_OPERATOR_IMPL

        std::string operator+(paramproxy const & p, char const * s) {
            return p.cast<std::string>() + s;
        }

        std::string operator+(char const * s, paramproxy const & p) {
            return s + p.cast<std::string>();
        }

    }
}
