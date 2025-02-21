/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
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

#ifndef OSIRIS_STD_IMPL_H
#define OSIRIS_STD_IMPL_H

// #include <palm/config.h>
#include <alps/osiris/dump.h>

#include <algorithm>
#include <complex>

namespace alps {
namespace detail {

/// function object for serialization
template<class T>
class saveTo {
  alps::ODump& dump_;
public:
  /// the constructor takes the ODump as argument
  saveTo(alps::ODump& dump) : dump_(dump) {}
  /// the function call operator serializes the object
  inline void operator()(const T& x) { dump_ << x;}
};
  
/// serialize a container
template<class C>
inline void saveContainer(alps::ODump& dump, const C& x)
{
  saveTo<typename C::value_type> save(dump);
  dump << uint32_t(x.size());
  std::for_each(x.begin(),x.end(),save);
}

/// deserialize a set-like container
template<class C>
inline void loadContainer(alps::IDump& dump, C& x)
{
  x=C();
  uint32_t n(dump);
  typename C::value_type val;
  while (n--) {
    dump >> val;
    x.push_back(val);
  }
}

/// deserialize a list-like container
template<class C>
inline alps::IDump& loadArrayLikeContainer(alps::IDump& dump, C& x)
{
  x.resize(uint32_t(dump));
  for(typename C::iterator p=x.begin();p!=x.end();++p)
    dump >> *p;
  return dump;
}

/// deserialize a stack-like container
template<class C>
inline alps::IDump& loadStackLikeContainer(alps::IDump& dump, C& x)
{
  x=C(); // empty stack
  int n(dump); // number of elements
  
  typename C::value_type elem;
  while (n--)
    {
      dump >> elem;
      x.push(elem);
    }
  return dump;
}

/// deserialize a stack-like container
template<class C>
inline alps::IDump& loadSetLikeContainer(alps::IDump& dump, C& x)
{
  x=C(); // empty stack
  int n(dump); // number of elements

  typename C::value_type elem;
  while (n--)
    {
      dump >> elem;
      x.insert(x.end(),elem);
    }
  return dump;
}


/** traits class specifying for whoch types there are optimized
    saveArray functions in the dump classes */

template <class T> struct TypeDumpTraits {
  BOOST_STATIC_CONSTANT(bool, hasArrayFunction=false);
};

#define _HAS_ARRAY(T) template<> struct TypeDumpTraits< T > {\
  BOOST_STATIC_CONSTANT(bool, hasArrayFunction=true);};
  
#ifndef BOOST_NO_INT64_T
_HAS_ARRAY(int64_t)
_HAS_ARRAY(uint64_t)
#endif
_HAS_ARRAY(int32_t)
_HAS_ARRAY(int16_t)
_HAS_ARRAY(int8_t)
_HAS_ARRAY(uint32_t)
_HAS_ARRAY(uint16_t)
_HAS_ARRAY(uint8_t)
_HAS_ARRAY(float)
_HAS_ARRAY(double)
_HAS_ARRAY(long double)
_HAS_ARRAY(bool)
_HAS_ARRAY(std::complex<float>)
_HAS_ARRAY(std::complex<double>)
_HAS_ARRAY(std::complex<long double>)

#undef _HAS_ARRAY

} // end namespace detail
} // end namespace alps

#endif // OSIRIS_std_IMPL_H
