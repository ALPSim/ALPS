/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
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
