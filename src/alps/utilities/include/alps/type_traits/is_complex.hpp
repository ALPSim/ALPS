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

#ifndef ALPS_TYPE_TRAITS_IS_COMPLEX_H
#define ALPS_TYPE_TRAITS_IS_COMPLEX_H

#include <boost/mpl/bool.hpp>
#include <complex>

// maybe we can automate this by checking for the existence of a value_type member

namespace alps {

template <class T>
struct is_complex : public boost::mpl::false_ {};

template <class T>
struct is_complex<std::complex<T> > : public boost::mpl::true_ {};

} // end namespace alps

#endif // ALPS_TYPE_TRAITS_IS_COMPLEX_H
