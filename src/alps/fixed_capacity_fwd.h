/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2003 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_FIXED_CAPACITY_FWD_H
#define ALPS_FIXED_CAPACITY_FWD_H

#include <cstddef>

namespace alps {

namespace fixed_capacity {

// Forward declarations of checking policy classes
//   definitions are given in <alps/fixed_capacity/checking.h>

struct no_checking;
struct capacity_checking;
struct strict_checking;

} // namespace fixed_capacity

// Forward declarations of alps::fixed_capacity_[vector,deque]
//   definitions are given in <alps/fixed_capacity_vector.h> and
//   <alps/fixed_capacity_deque.h>, respectively

template<class T,
         std::size_t N,
         class CheckingPolicy = ::alps::fixed_capacity::no_checking>
class fixed_capacity_vector;
template<class T,
         std::size_t N,
         class CheckingPolicy = ::alps::fixed_capacity::no_checking>
class fixed_capacity_deque;

// Forward declaration of traits class alps::fixed_capacity_traits
//   definition is given in <alps/fixed_capacity_traits.h>

template<class C> struct fixed_capacity_traits;

} // namespace alps

#endif // ALPS_FIXED_CAPACITY_FWD_H
