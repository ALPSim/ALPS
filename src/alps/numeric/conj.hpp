/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2012 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>,
*                            Andreas Hehn <hehn@phys.ethz.ch>                   
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

/* $Id$ */

#ifndef ALPS_NUMERIC_CONJ_HPP
#define ALPS_NUMERIC_CONJ_HPP

#include <alps/numeric/matrix/entity.hpp>
#include <complex>
#include <boost/type_traits/is_fundamental.hpp>
#include <boost/static_assert.hpp>

namespace alps { namespace numeric {


template <class T>
typename boost::enable_if<boost::is_fundamental<T>,T>::type conj (T x)
{ 
  return x;
}

// if std::complex<T> is used std::conj will be called by argument dependent look-up

template <class T>
typename boost::enable_if<boost::is_fundamental<T>,void>::type conj_inplace(T& t, tag::scalar)
{
}

template <class T>
void conj_inplace(std::complex<T>& x, tag::scalar)
{
  BOOST_STATIC_ASSERT((boost::is_fundamental<T>::value));
  using std::conj;
  x = conj(x);
}

template <typename T>
void conj_inplace(T& t)
{
    conj_inplace(t, typename get_entity<T>::type());
}

} }  // end namespace alps::numeric

#endif // ALPS_NUMERIC_CONJ_HPP
