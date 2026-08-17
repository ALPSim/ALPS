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
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
