/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Matthias Troyer <troyer@comp-phys.org>,
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
