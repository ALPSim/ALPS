/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2012        by Michele Dolfi <dolfim@phys.ethz.ch>                *
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

#include <iostream>
#include <iterator>
#include <complex>
#include <vector>
#include <alps/numeric/real.hpp>

#define BOOST_TEST_SOURCE
#define BOOST_TEST_MODULE alps::numeric::real
#ifndef ALPS_LINK_BOOST_TEST
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif
#include <boost/mpl/list.hpp>

//
// List of types T for which the real(T) is tested
//
typedef boost::mpl::list<
      float
    , double
    , std::complex<float>
    , std::complex<double>
    , std::vector<float>
    , std::vector<double>
    , std::vector<std::complex<float> >
    , std::vector<std::complex<double> >
    , std::vector<std::vector<float> >
    , std::vector<std::vector<double> >
    , std::vector<std::vector<std::complex<float> > >
    , std::vector<std::vector<std::complex<double> > >
> test_types;


//
// ostream overloads
//
// template <typename T>
// std::ostream& operator<< (std::ostream& os, std::vector<T const&> vec)
// {
//     os << "[";
//     std::copy(vec.begin(), vec.end(), std::ostream_iterator<T const&>(os," "));
//     os << "]";
//     return os;
// }
template <typename T>
std::ostream& operator<< (std::ostream& os, std::vector<T> const & vec)
{
    os << "[";
    for (std::size_t i=0; i<vec.size(); ++i)
        os << (i ? " " : "") << vec[i];
    os << "]";
    return os;
}


//
// Filling functions
//
float fill_val = 1;
const std::size_t vecsize = 2;
template <typename T>
void fill (T & v)
{
    v = (fill_val++);
}
template <typename T>
void fill (std::complex<T> & v) {
    T real_part = fill_val++;
    T imag_part = fill_val++;
    v = std::complex<T>(real_part,imag_part);
}
template <typename T>
void fill (std::vector<T> & v)
{
    v.resize(vecsize);
    std::for_each(v.begin(), v.end(), static_cast<void (*)(T &)>(&fill));
}



//
// Test with full namespaces
//
BOOST_AUTO_TEST_CASE_TEMPLATE( real_with_namespace, T, test_types )
{
    T val; fill(val);
    std::cout << "real( " << val << " ) = " << alps::numeric::real(val) << std::endl;
}

//
// Test letting the compiler resolve the overloads
//
BOOST_AUTO_TEST_CASE_TEMPLATE( real_without_namespace, T, test_types )
{
    T val; fill(val);
    using alps::numeric::real;
    std::cout << "real( " << val << " ) = " << real(val) << std::endl;
}

