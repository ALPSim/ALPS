/***************************************************************************
* ALPS library
*
* fixed_capacity/checking.h header file
*
* $Id$
*
* Copyright (C) 2002-2003 by Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person or organization 
* obtaining a copy of the software covered by this license (the "Software") 
* to use, reproduce, display, distribute, execute, and transmit the Software, 
* and to prepare derivative works of the Software, and to permit others
* to do so for non-commerical academic use, all subject to the following:
*
* The copyright notice in the Software and this entire statement, including 
* the above license grant, this restriction and the following disclaimer, 
* must be included in all copies of the Software, in whole or in part, and 
* all derivative works of the Software, unless such copies or derivative 
* works are solely in the form of machine-executable object code generated by 
* a source language processor.

* In any scientific publication based in part or wholly on the Software, the
* use of the Software has to be acknowledged and the publications quoted
* on the web page http://www.alps.org/license/ have to be referenced.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
**************************************************************************/

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
