/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Lukas Gamper <gamperl@gmail.com>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id: funcitonal.hpp 3958 2010-03-05 09:24:06Z gamperl $ */

#ifndef ALPS_NUMERIC_FUNCTIONAL_HPP
#define ALPS_NUMERIC_FUNCTIONAL_HPP

#include <alps/numeric/vector_functions.hpp>
#include <alps/boost/accumulators/numeric/functional.hpp>
#include <alps/boost/accumulators/numeric/functional/vector.hpp>

namespace alps { 
    namespace numeric {
        template <typename T> struct unary_minus : public std::unary_function<T, T> {
            T operator()(T const & x) const {
                using boost::numeric::operators::operator-;
                return -x; 
            }
        };

        template <typename T, typename U, typename R> struct plus : public std::binary_function<T, U, R> {
            R operator()(T const & x, U const & y) const {
                using boost::numeric::operators::operator+;
                return x + y; 
            }
        };
        template <typename T> struct plus<T, T, T> : public std::binary_function<T, T, T> {
            T operator()(T const & x, T const & y) const {
                using boost::numeric::operators::operator+;
                return x + y; 
            }
        };

        template <typename T, typename U, typename R> struct minus : public std::binary_function<T, U, R> {
            R operator()(T const & x, U const & y) const {
                using boost::numeric::operators::operator-;
                return x - y; 
            }
        };
        template <typename T> struct minus<T, T, T> : public std::binary_function<T, T, T> {
            T operator()(T const & x, T const & y) const {
                using boost::numeric::operators::operator-;
                return x - y; 
            }
        };

        template <typename T, typename U, typename R> struct multiplies : public std::binary_function<T, U, R> {
            R operator()(T const & x, U const & y) const {
                using boost::numeric::operators::operator*;
                return x * y; 
            }
        };
        template <typename T> struct multiplies<T, T, T> : public std::binary_function<T, T, T> {
            T operator()(T const & x, T const & y) const {
                using boost::numeric::operators::operator*;
                return x * y; 
            }
        };

        template <typename T, typename U, typename R> struct divides : public std::binary_function<T, U, R> {
            R operator()(T const & x, U const & y) const {
                using boost::numeric::operators::operator/;
                return x / y; 
            }
        };
        template <typename T> struct divides<T, T, T> : public std::binary_function<T, T, T> {
            T operator()(T const & x, T const & y) const {
                using boost::numeric::operators::operator/;
                return x / y; 
            }
        };
    } 
}

#endif // ALPS_NUMERIC_FUNCTIONAL_HPP
