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

#ifndef ALPS_TYPE_TRAITS_IS_SEQUENCE_H
#define ALPS_TYPE_TRAITS_IS_SEQUENCE_H

#include <alps/config.h>
#include <boost/mpl/bool.hpp>
#include <alps/type_traits/has_value_type.hpp>
#include <valarray>
#include <vector>
#include <complex>

// maybe we can automate this by checking for the existence of a value_type member

namespace alps {

template <class T>
struct is_sequence : public alps::has_value_type<T> {};

template <class T>
struct is_sequence<std::valarray<T> > : public boost::mpl::true_ {};

template <class T>
struct is_sequence<std::complex<T> > : public boost::mpl::false_ {};
 
template <>
struct is_sequence<std::string> : public boost::mpl::false_ {};

} // end namespace alps

#endif // ALPS_TYPE_TRAITS_ELEMENT_TYPE_H
