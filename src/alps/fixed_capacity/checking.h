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

#ifndef ALPS_FIXED_CAPACITY_CHECKING_H
#define ALPS_FIXED_CAPACITY_CHECKING_H

#include <cstddef>
#include <stdexcept>

namespace alps {
  
namespace fixed_capacity {

// struct no_checking -------------------------------------------------------//

struct no_checking
{
  // for fixed_capacity_vector
  template<std::size_t N>
  struct vector {
    static void range_check(std::size_t, std::size_t) {}
    static void capacity_check(std::size_t) {}
    template<class T>
    static void pointer_check(const T*, const T*, const T*) {}
  };

  // for fixed_capacity_deque
  template<std::size_t N>
  struct deque {
    static void range_check(std::size_t, std::size_t) {}
    static void capacity_check(std::size_t) {}
    template<class T>
    static void pointer_check(const T*, const T*, const T*, const T*) {}
  };

}; // no_checking


// struct capacity_checking -------------------------------------------------//

struct capacity_checking
{
  // for fixed_capacity_vector
  template<std::size_t N>
  struct vector {
    static void range_check(std::size_t, std::size_t) {}
    static void capacity_check(std::size_t n) {
      if (n > N) throw std::range_error("fixed_capacity_vector");
    }
    template<class T>
    static void pointer_check(const T*, const T*, const T*) {}
  };

  // for fixed_capacity_deque
  template<std::size_t N>
  struct deque {
    static void range_check(std::size_t, std::size_t) {}
    static void capacity_check(std::size_t n) {
      if (n > N) throw std::range_error("fixed_capacity_deque");
    }
    template<class T>
    static void pointer_check(const T*, const T*, const T*, const T*) {}
  };

}; // capacity_checking


// struct strict_checking ---------------------------------------------------//

struct strict_checking
{
  // for fixed_capacity_vector
  template<std::size_t N>
  struct vector {
    static void range_check(std::size_t s, std::size_t i) {
      if (i >= s) throw std::range_error("fixed_capacity_vector");
    }
    static void capacity_check(std::size_t n) {
      if (n > N) throw std::range_error("fixed_capacity_vector");
    }
    template<class T>
    static void pointer_check(const T* base, const T* last, const T* ptr) {
      if (ptr < base || ptr > last)
        throw std::range_error("fixed_capacity_vector");
    }
  };

  // for fixed_capacity_deque
  template<std::size_t N>
  struct deque {
    static void range_check(std::size_t s, std::size_t i) {
      if (i >= s) throw std::range_error("fixed_capacity_deque");
    }
    static void capacity_check(std::size_t n) {
      if (n > N) throw std::range_error("fixed_capacity_deque");
    }
    template<class T>
    static void pointer_check(const T* base, const T* first, const T* last,
                              const T* ptr) {
      if (last - first >= 0) {
        if (ptr < first || ptr > last) {
          throw std::range_error("fixed_capacity_deque");
        }
      } else {
        if (ptr < base || (ptr > last && ptr < first) || ptr > base + N) {
          throw std::range_error("fixed_capacity_deque");
        }
      }
    }
  };

}; // strict_checking

} // namespace fixed_capacity

} // namespace alps

#endif // ALPS_FIXED_CAPACITY_CHECKING_H
