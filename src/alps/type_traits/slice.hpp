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

#ifndef ALPS_TYPE_TRAITS_SLICE_HPP
#define ALPS_TYPE_TRAITS_SLICE_HPP

#include <alps/type_traits/is_sequence.hpp>
#include <alps/type_traits/element_type.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/lexical_cast.hpp>

namespace alps {

template <class T>
struct slice_index
{
  typedef std::size_t type;
};

template <class T>
typename boost::enable_if<is_sequence<T>,
  std::pair<typename slice_index<T>::type,typename slice_index<T>::type >
>::type
slices (T const& x) 
{ 
  return std::pair<typename slice_index<T>::type,typename slice_index<T>::type >(0,x.size());
}

template <class T>
typename boost::disable_if<is_sequence<T>,
  std::pair<typename slice_index<T>::type,typename slice_index<T>::type >
>::type
slices (T const&) 
{ 
  return std::pair<typename slice_index<T>::type,typename slice_index<T>::type >(0,1);
}

template <class ValueType>
typename boost::enable_if<is_sequence<ValueType>, std::string>::type
slice_name(ValueType const& ,unsigned i) 
{ 
  return boost::lexical_cast<std::string,int>(i);
}

template <class ValueType>
typename boost::disable_if<is_sequence<ValueType>, std::string>::type
slice_name(ValueType const& ,unsigned) 
{ 
  return "";
}

template <class ValueType>
typename boost::enable_if<is_sequence<ValueType>, 
  typename element_type<ValueType>::type>::type
slice_value(ValueType const& x ,unsigned i) 
{ 
  return i < x.size() ? x[i] : typename element_type<ValueType>::type();
}


template <class ValueType>
typename boost::disable_if<is_sequence<ValueType>, ValueType const&>::type
slice_value(ValueType const& x ,unsigned) 
{ 
  return x;
}

template <class ValueType>
typename boost::enable_if<is_sequence<ValueType>, 
  typename element_type<ValueType>::type&>::type
slice_value(ValueType& x,unsigned i) 
{ 
  return x[i];
}


template <class ValueType>
typename boost::disable_if<is_sequence<ValueType>, ValueType&>::type
slice_value(ValueType& x ,unsigned) 
{ 
  return x;
}

template <class ValueType>
struct slice_it  
{                      
  typedef typename element_type<ValueType>::type result_type;
  typedef ValueType const& first_argument_type;
  typedef typename slice_index<ValueType>::type second_argument_type;
   
  result_type operator()(ValueType const& x, second_argument_type i) const
  { 
    return slice_value(x,i);
  }
};                                             


} // end namespace alps

#endif // ALPS_TYPE_TRAITS_ELEMENT_TYPE_H
