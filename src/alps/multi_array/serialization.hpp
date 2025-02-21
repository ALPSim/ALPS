/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2012 by Ilia Zintchenko <iliazin@gmail.com>                       *
 *                       Jan Gukelberger                                           *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_MULTI_ARRAY_SERIALIZATION_HPP
#define ALPS_MULTI_ARRAY_SERIALIZATION_HPP

#include <alps/multi_array/multi_array.hpp>

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/serialization.hpp>

namespace alps{

  template<typename Archive, typename T, std::size_t D, class Allocator> 
  inline void save(Archive & ar, const multi_array<T,D,Allocator> & t, const unsigned int file_version) 
  { 
    using boost::serialization::make_nvp;
    using boost::serialization::make_array;
    ar << make_nvp("dimensions", make_array(t.shape(), D)); 
    ar << make_nvp("data", make_array(t.data(), t.num_elements()));
  } 

  template<typename Archive, typename T, std::size_t D, class Allocator> 
  inline void load(Archive & ar, multi_array<T,D,Allocator> & t, const unsigned int file_version) 
  { 
    using boost::serialization::make_nvp;
    using boost::serialization::make_array;

    typedef typename multi_array<T,D,Allocator>::size_type size_type;

    boost::array<size_type, D> dimensions;
    ar >> make_nvp("dimensions", make_array(dimensions.c_array(), D));
    t.resize(dimensions); 
    ar >> make_nvp("data", make_array(t.data(), t.num_elements())); 
  } 

  template<typename Archive, typename T, std::size_t D, class Allocator> 
  inline void serialize(Archive & ar, multi_array<T,D,Allocator>& t, const unsigned int file_version) 
  { 
    using boost::serialization::split_free;
    
    split_free(ar, t, file_version); 
  }

}//namespace alps

#endif // ALPS_MULTI_ARRAY_SERIALIZATION_HPP
