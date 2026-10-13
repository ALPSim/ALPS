/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_TYPE_TRAITS_ITERATOR_TYPE_HPP
#define ALPS_TYPE_TRAITS_ITERATOR_TYPE_HPP

#include <valarray>

namespace alps {

template <class Collection> 
struct iterator_type
{
  typedef typename Collection::iterator type;
};

template <class Collection> 
struct const_iterator_type
{
  typedef typename Collection::iterator type;
};

template <class T> 
struct iterator_type<std::valarray<T> >
{
  typedef T* type;
};

template <class T> 
struct const_iterator_type<std::valarray<T> >
{
  typedef T const * type;
};

} // namespace alps

#endif // ALPS_TYPE_TRAITS_ITERATOR_TYPE_HPP
