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

#ifndef ALPS_TYPE_TRAITS_ELEMENT_TYPE_H
#define ALPS_TYPE_TRAITS_ELEMENT_TYPE_H

#include <alps/type_traits/has_value_type.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/bool.hpp>
 
 namespace alps {
 
template <class T> struct element_type_recursive;

namespace detail {

template <class T, class F>
struct element_type_helper {};
  
 template <class T>
struct element_type_helper<T,boost::mpl::false_> 
{
  typedef T type;
};

template <class T>
struct element_type_helper<T,boost::mpl::true_> 
{
    typedef typename T::value_type type;
};

template <class T, class F>
struct element_type_recursive_helper {};
  
 template <class T>
struct element_type_recursive_helper<T,boost::mpl::false_> 
{
  typedef T type;
};

template <class T>
struct element_type_recursive_helper<T,boost::mpl::true_> 
    : element_type_recursive<typename T::value_type>
{
};


}

template <class T>
 struct element_type
 : public detail::element_type_helper<T,typename has_value_type<T>::type > {};

 template <class T>
 struct element_type_recursive
 : public detail::element_type_recursive_helper<T,typename has_value_type<T>::type > {};

}
 
#endif // ALPS_TYPE_TRAITS_ELEMENT_TYPE_H
