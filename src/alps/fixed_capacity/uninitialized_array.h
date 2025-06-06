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

#ifndef ALPS_FIXED_CAPACITY_UNINITIALIZED_ARRAY_HPP
#define ALPS_FIXED_CAPACITY_UNINITIALIZED_ARRAY_HPP

#include <boost/config.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/alignment_traits.hpp>
#include <cstddef>

namespace alps {

// class template uninitialized_array ---------------------------------------//

template<class T, std::size_t N>
class uninitialized_array
{
public:
  // types:
  typedef std::size_t size_type;
  typedef T           value_type;
  typedef T&          reference;
  typedef const T&    const_reference;
  typedef T*          iterator;
  typedef const T*    const_iterator;

  BOOST_STATIC_CONSTANT(size_type, static_size = N);

  // compiler-generated constructors/destructor are fine

  // iterators:
  iterator begin() { return reinterpret_cast<iterator>(buffer_); }
  const_iterator begin() const {
    return reinterpret_cast<const_iterator>(buffer_);
  }
  iterator end() { return begin() + N; }
  const_iterator end() const { return begin() + N; }

  // capacity:
  static size_type size() { return N; }

  // element access:
  reference operator[](size_type i) { return *(begin() + i); }
  const_reference operator[](size_type i) const { return *(begin() + i); }
  
private:
  BOOST_STATIC_ASSERT(N > 0);

  union {
    char buffer_[N * sizeof(T)];
    typename boost::type_with_alignment<boost::alignment_of<T>::value>::type
      dummy_;
  };

}; // uninitialized_array

} // end namespace alps

#endif // ALPS_FIXED_CAPACITY_UNINITIALIZED_ARRAY_HPP
