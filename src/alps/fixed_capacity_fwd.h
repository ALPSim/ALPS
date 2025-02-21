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
