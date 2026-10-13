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

#ifndef ALPS_TYPE_TRAITS_COVARIANCE_TYPE_H
#define ALPS_TYPE_TRAITS_COVARIANCE_TYPE_H

#include <alps/type_traits/average_type.hpp>
#include <alps/type_traits/element_type.hpp>
#include <alps/type_traits/is_sequence.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/mpl/if.hpp>


namespace alps {

template <class T>
struct covariance_type
{
 typedef typename boost::mpl::if_<
     is_sequence<T>,
     typename boost::numeric::ublas::matrix<
       typename average_type<typename element_type<T>::type>::type
     >,
     typename average_type<T>::type
   >::type type;
};


} // end namespace alps

#endif // ALPS_TYPE_TRAITS_COVARIANCE_TYPE_H
