/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Lukas Gamper <gamperl@gmail.com>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
        template <typename T> struct unary_minus {
            T operator()(T const & x) const {
                using boost::numeric::operators::operator-;
                return -x;
            }
        };

        template <typename T, typename U, typename R> struct plus {
            R operator()(T const & x, U const & y) const {
                using boost::numeric::operators::operator+;
                return x + y;
            }
        };
        template <typename T> struct plus<T, T, T> {
            T operator()(T const & x, T const & y) const {
                using boost::numeric::operators::operator+;
                return x + y;
            }
        };

        template <typename T, typename U, typename R> struct minus {
            R operator()(T const & x, U const & y) const {
                using boost::numeric::operators::operator-;
                return x - y;
            }
        };
        template <typename T> struct minus<T, T, T> {
            T operator()(T const & x, T const & y) const {
                using boost::numeric::operators::operator-;
                return x - y;
            }
        };

        template <typename T, typename U, typename R> struct multiplies {
            R operator()(T const & x, U const & y) const {
                using boost::numeric::operators::operator*;
                return x * y;
            }
        };
        template <typename T> struct multiplies<T, T, T> {
            T operator()(T const & x, T const & y) const {
                using boost::numeric::operators::operator*;
                return x * y;
            }
        };

        template <typename T, typename U, typename R> struct divides {
            R operator()(T const & x, U const & y) const {
                using boost::numeric::operators::operator/;
                return x / y;
            }
        };
        template <typename T> struct divides<T, T, T> {
            T operator()(T const & x, T const & y) const {
                using boost::numeric::operators::operator/;
                return x / y;
            }
        };
    } 
}

#endif // ALPS_NUMERIC_FUNCTIONAL_HPP
