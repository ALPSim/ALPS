/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2003 by Synge Todo <wistaria@comp-phys.org>
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

/* $Id$ */

#ifndef ALPS_FIXED_CAPACITY_TRAITS_H
#define ALPS_FIXED_CAPACITY_TRAITS_H

#include <cstddef>
#include <queue>
#include <stack>
#include <boost/config.hpp>
#include <alps/fixed_capacity_fwd.h>

namespace alps {

// traits class alps::fixed_capacity_traits --------------------------------//

template<class C>
struct fixed_capacity_traits {
  BOOST_STATIC_CONSTANT(bool, capacity_is_fixed = false);
};

// specializations for fixed_capacity_[vector,deque]

template<class T, std::size_t N, class C>
struct fixed_capacity_traits<fixed_capacity_vector<T, N, C> > {
  BOOST_STATIC_CONSTANT(bool, capacity_is_fixed = true);
  BOOST_STATIC_CONSTANT(std::size_t, static_max_size = N);
};


template<class T, std::size_t N, class C>
struct fixed_capacity_traits<fixed_capacity_deque<T, N, C> > {
  BOOST_STATIC_CONSTANT(bool, capacity_is_fixed = true);
  BOOST_STATIC_CONSTANT(std::size_t, static_max_size = N);
};

// specializations for adaptors using fixed_capacity_[vector,deque] as
// a base container

template<class T, class C>
struct fixed_capacity_traits<std::stack<T, C> >
  : public fixed_capacity_traits<C> {};

template<class T, class C>
struct fixed_capacity_traits<std::queue<T, C> >
  : public fixed_capacity_traits<C> {};

template<class T, class C, class Cmp>
struct fixed_capacity_traits<std::priority_queue<T, C, Cmp> >
  : public fixed_capacity_traits<C> {};

#ifndef BOOST_NO_INCLASS_MEMBER_INITIALIZATION
template <class C>
const bool fixed_capacity_traits<C>::capacity_is_fixed;

template<class T, std::size_t N, class C>
const bool fixed_capacity_traits<fixed_capacity_vector<T, N, C> >::capacity_is_fixed;

template<class T, std::size_t N, class C>
const std::size_t fixed_capacity_traits<fixed_capacity_vector<T, N, C> >::static_max_size;

template<class T, std::size_t N, class C>
const bool fixed_capacity_traits<fixed_capacity_deque<T, N, C> >::capacity_is_fixed;

template<class T, std::size_t N, class C>
const std::size_t fixed_capacity_traits<fixed_capacity_deque<T, N, C> >::static_max_size;
#endif

} // namespace alps

#endif // ALPS_FIXED_CAPACITY_TRAITS_H
