/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010        by Andreas Hehn <hehn@phys.ethz.ch>                   *
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

#ifndef MATRIX_UNIT_TESTS_HPP
#define MATRIX_UNIT_TESTS_HPP

#define BOOST_TEST_SOURCE
#define BOOST_TEST_MODULE alps::numeric::matrix

#ifndef ALPS_LINK_BOOST_TEST
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif
#include <boost/filesystem.hpp>
#include <boost/mpl/list.hpp>

#include <boost/lambda/lambda.hpp>
#include <complex>
#include <numeric>
#include <iostream>

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/matrix_interface.hpp>
#include <alps/numeric/matrix/vector.hpp>

//
// List of types T for which the matrix<T> is tested
//
typedef boost::mpl::list<
      float
    , double
    , int
    , unsigned int
    , long unsigned int
    , std::complex<float>
    , std::complex<double>
> test_types;
// long long unsigned int causes problems in boost::iterator facade


namespace type_pairs
{
struct DComplexDouble
{
    typedef std::complex<double> first_type;
    typedef double second_type;
    typedef std::complex<double> result_type;
};

struct DoubleDComplex
{
    typedef double first_type;
    typedef std::complex<double> second_type;
    typedef std::complex<double> result_type;
};

struct IntDouble
{
    typedef int first_type;
    typedef double second_type;
    typedef double result_type;
};

struct DoubleInt
{
    typedef double first_type;
    typedef int second_type;
    typedef double result_type;
};
};

//
// List of type pairs <T,U> for which the mixed type matrix std::vector multiplication is tested.
//
typedef boost::mpl::list<type_pairs::IntDouble, type_pairs::DoubleInt, type_pairs::DoubleDComplex, type_pairs::DComplexDouble> test_type_pairs;

template <typename OutputIterator, typename T>
T fill_range_with_numbers(OutputIterator begin, OutputIterator end, T iota)
{
    // Unfortunately we can't use the postincrement operator, due to std:complex<>
    // -> so we have to emulate it's behaviour...
    std::transform(begin,end,begin,boost::lambda::_1 = (boost::lambda::var(iota)+=T(1))-T(1));
    return iota;
}

template <typename T>
T fill_matrix_with_numbers(alps::numeric::matrix<T>& a, T iota = T(0))
{
    for(unsigned int i=0; i<num_rows(a); ++i)
    {
        std::pair<typename alps::numeric::matrix<T>::row_element_iterator, typename alps::numeric::matrix<T>::row_element_iterator> range(row(a,i));
        iota += fill_range_with_numbers(range.first,range.second,T(i));
    }
    return iota;
}

#endif // MATRIX_UNIT_TESTS_HPP
