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

#ifndef ALPS_NGS_DETAIL_GET_NUMPY_TYPE_HPP
#define ALPS_NGS_DETAIL_GET_NUMPY_TYPE_HPP

#include <alps/ngs/config.hpp>

#if !defined(ALPS_HAVE_PYTHON)
    #error numpy is only available if python is enabled
#endif

#include <alps/ngs/boost_python.hpp>
#include <alps/python/numpy_import.hpp>

#include <complex>

#define ALPS_NGS_FOREACH_NATIVE_NUMPY_TYPE(CALLBACK)                                                                                            \
    CALLBACK(bool)                                                                                                                              \
    CALLBACK(char)                                                                                                                              \
    CALLBACK(signed char)                                                                                                                       \
    CALLBACK(unsigned char)                                                                                                                     \
    CALLBACK(short)                                                                                                                             \
    CALLBACK(unsigned short)                                                                                                                    \
    CALLBACK(int)                                                                                                                               \
    CALLBACK(unsigned)                                                                                                                          \
    CALLBACK(long)                                                                                                                              \
    CALLBACK(unsigned long)                                                                                                                     \
    CALLBACK(long long)                                                                                                                         \
    CALLBACK(unsigned long long)                                                                                                                \
    CALLBACK(float)                                                                                                                             \
    CALLBACK(double)                                                                                                                            \
    CALLBACK(long double)                                                                                                                       \
    CALLBACK(std::complex<float>)                                                                                                               \
    CALLBACK(std::complex<double>)                                                                                                              \
    CALLBACK(std::complex<long double>)

namespace alps {
    namespace detail {
        #define ALPS_NGS_DECL_NUMPY_TYPE(T)                                                                                                     \
            ALPS_DECL int get_numpy_type(T);
        ALPS_NGS_FOREACH_NATIVE_NUMPY_TYPE(ALPS_NGS_DECL_NUMPY_TYPE)
        #undef ALPS_NGS_DECL_NUMPY_TYPE
    }
}

#endif
