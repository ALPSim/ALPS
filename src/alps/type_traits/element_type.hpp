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
