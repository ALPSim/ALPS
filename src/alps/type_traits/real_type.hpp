/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Matthias Troyer <troyer@comp-phys.org>,
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

/* $Id: obsvalue.h 3435 2009-11-28 14:45:38Z troyer $ */

#ifndef ALPS_TYPE_TRAITS_REAL_TYPE_H
#define ALPS_TYPE_TRAITS_REAL_TYPE_H

#include <boost/mpl/bool.hpp>
#include <complex>
#include <valarray>
#include <vector>

// maybe we can automate this by checking for the existence of a value_type member

namespace alps {

template <class T>
struct real_type 
{
  typedef T type;
};

template <class T>
struct real_type<std::complex<T> > : public real_type<T> {};


template <class T>
struct real_type<std::valarray<T> > {
  typedef std::valarray<typename real_type<T>::type> type;
};

template <class T, class A>
struct real_type<std::vector<T,A> > {
  typedef std::vector<typename real_type<T>::type,A> type;
};



} // end namespace alps

#endif // ALPS_TYPE_TRAITS_REAL_TYPE_H
