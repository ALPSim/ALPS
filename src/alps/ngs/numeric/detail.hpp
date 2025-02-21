/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Mario Koenz <mkoenz@ethz.ch>                       *
 * Copyright (C) 2012 - 2014 by Lukas Gamper <gamperl@gmail.com>                   *
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

#ifndef ALPS_NGS_NUMERIC_DETAIL_HEADER
#define ALPS_NGS_NUMERIC_DETAIL_HEADER

#include <alps/ngs/stacktrace.hpp>
#include <alps/utility/resize.hpp>

#include <alps/multi_array.hpp>

#include <boost/array.hpp>

#include <vector>
#include <stdexcept>

namespace alps {
    namespace ngs { //merged with alps/numerics/vector_function.hpp
        namespace numeric {
            namespace detail {

                template<typename T, typename U>
                inline void check_size(T & a, U const & b) {}

                template<typename T, typename U>
                inline void check_size(std::vector<T> & a, std::vector<U> const & b) {
                    if(a.size() == 0)
                        alps::resize_same_as(a, b);
                    else if(a.size() != b.size())
                        boost::throw_exception(std::runtime_error("vectors must have the same size!" + ALPS_STACKTRACE));
                }

                template<typename T, typename U, std::size_t N, std::size_t M>
                inline void check_size(boost::array<T, N> & a, boost::array<U, M> const & b) {
                    boost::throw_exception(std::runtime_error("boost::array s must have the same size!" + ALPS_STACKTRACE));
                }

                template<typename T, typename U, std::size_t N>
                inline void check_size(boost::array<T, N> & a, boost::array<U, N> const & b) {}

                template<typename T, typename U, std::size_t D>
                inline void check_size(alps::multi_array<T, D> & a, alps::multi_array<U, D> const & b) {}
                
            }
        }
    }
}

#endif
