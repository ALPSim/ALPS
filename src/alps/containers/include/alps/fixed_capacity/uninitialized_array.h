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
