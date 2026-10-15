// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/vector.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#ifndef __DMTK_VECTOR_H__
#define __DMTK_VECTOR_H__

#include "dmtk/output.hpp"
#include <iostream>
#include <vector>
#include <iosfwd>
#include "dmtk/conj.hpp"
#include "dmtk/slice_iter.hpp"
#include "dmtk/gslice_iter.hpp"
#include "dmtk/range.hpp"
#include "dmtk/array_util.hpp"

namespace dmtk
{

// Vector<T> is a thin std::vector<T> wrapper that adds:
//   - dmtk-style \c operator()(size_t) and \c operator()(Range) accessors
//   - element-wise arithmetic operators
//   - binary read()/write() to std::istream/ostream
//   - randomize() seeded for ASLR-independent reproducibility
//
// Earlier versions kept a parallel quartet of \c v_start/v_end/v_size/v_capacity
// raw fields synchronised by an init() helper. That bookkeeping is redundant
// once the underlying container is std::vector — std::vector::data()/size()
// give us the same information for free.

template<class T>
class Vector : public std::vector<T>
{
  public:
    using vector = std::vector<T>;
    using iterator = T*;
    using const_iterator = const T*;

    Vector() = default;

    explicit Vector(size_t s) : vector(s) {}
    Vector(size_t s, const T& v) : vector(s, v) {}
    Vector(size_t s, const T* v) : vector(s, v) {}

    Vector(const Vector<T>& v) : vector(v) {}

    Vector(cslice_iter<T> v) : vector(v.size())
      {
        const size_t n = std::min(this->size(), v.size());
        array_copy2(n, v, *this);
      }

    Vector(slice_iter<T> v) : vector(v.size())
      {
        const size_t n = std::min(this->size(), v.size());
        array_copy2(n, v, *this);
      }

//  Operators

    T const& operator()(size_t i) const { return vector::operator[](i); }
    T&       operator()(size_t i)       { return vector::operator[](i); }

    Vector& operator=(const Vector<T>&);
    Vector& operator=(slice_iter<T>);
    Vector& operator=(cslice_iter<T>);
    Vector& operator=(const_iterator);
    Vector& operator=(const T&);

    Vector& operator+=(const Vector<T> &);
    Vector& operator-=(const Vector<T> &);
    Vector& operator*=(const Vector<T> &);
    Vector& operator/=(const Vector<T> &);

    Vector& operator+=(const T &);
    Vector& operator-=(const T &);
    Vector& operator*=(const T &);
    Vector& operator/=(const T &);

//  Iterator

    T* array()                   { return vector::data(); }
    const T* array() const       { return vector::data(); }
    iterator begin()             { return vector::data(); }
    const_iterator begin() const { return vector::data(); }
    iterator end()               { return vector::data() + this->size(); }
    const_iterator end() const   { return vector::data() + this->size(); }

//  Ranges

    slice_iter<T> operator()(Range);
    cslice_iter<T> operator()(Range) const;

//  Memory management

    inline size_t memory() const { return this->capacity() * sizeof(T); }
    inline size_t usage()  const { return this->size() * sizeof(T); }
    Vector& resize(size_t n)            { vector::resize(n);    return *this; }
    Vector& resize(size_t n, const T& v){ vector::resize(n, v); return *this; }
    Vector& reserve(size_t n)           { vector::reserve(n);   return *this; }

//  Methods

    // The owned invocation saves, resets and restores this per-scalar stream.
    // Unowned CLI calls keep the original process-lifetime sequence.
    static long& random_state() noexcept
    {
      static long state = 42;
      return state;
    }

    Vector<T>& randomize()
    {
      // Fixed seed for ASLR-independent reproducibility. The previous
      // `abs((long)this)/0xffff` expression invoked pointer->long cast
      // UB and made Lanczos random-start vectors change per process run,
      // forcing loose (>= 1e-6) tolerances on dmrg numerical locks.
      long& idum = random_state();
      for (size_t i = 0; i < this->size(); i++)
        operator()(i) = quickran(idum);
      return *this;
    }

//  Streams

    void read(std::istream& s);
    void write(std::ostream& s) const;
};


template<class T>
slice_iter<T>
Vector<T>::operator()(Range range)
{
  std::slice s(range.start(), range.size(), range.stride());
  slice_iter<T> iter(this, s);
  return iter;
}


template<class T>
cslice_iter<T>
Vector<T>::operator()(Range range) const
{
  std::slice s(range.start(), range.size(), range.stride());
  cslice_iter<T> iter(this, s);
  return iter;
}


template<class T>
inline Vector<T>&
Vector<T>::operator=(const Vector<T> &v)
{
  vector::operator=(v);
  return *this;
}

template<class T>
inline Vector<T>&
Vector<T>::operator=(const T& v)
{
  std::fill(this->begin(), this->end(), v);
  return *this;
}

template<class T>
inline Vector<T>&
Vector<T>::operator=(const_iterator v)
{
  array_copy(this->size(), v, this->begin());
  return *this;
}

template<class T>
inline Vector<T>&
Vector<T>::operator=(slice_iter<T> v)
{
  resize(v.size());
  array_copy2(v.size(), v, *this);
  return *this;
}

template<class T>
inline Vector<T>&
Vector<T>::operator=(cslice_iter<T> v)
{
  resize(v.size());
  array_copy2(v.size(), v, *this);
  return *this;
}

////////////////////////////////////////////////////////

#define BINARY_OP(op,ap) \
template<class T> \
inline Vector<T>& \
Vector<T>::op(const Vector<T>& mm) \
{ \
  iterator dest = begin(); \
  const_iterator orig = mm.begin(); \
  size_t n = this->size(); \
  while(n--) *dest++ ap *orig++; \
  return *this; \
}

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline Vector<T>& \
Vector<T>::op(const T& v) \
{ \
  iterator pos = begin(); \
  while(pos != end()) *pos++ ap v; \
  return *this; \
}

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

////////////////////////////////////////////////////////
// Materialising arithmetic operators. Prior versions of this header
// returned boost.ublas-style expression templates via meta.h; the
// operators below allocate and return concrete Vector<T> instances
// instead. Subclasses (dmtk::VectorState) provide their own overloads
// returning the derived type and will shadow these via overload
// resolution when their operands are more specific.

template<class T>
inline Vector<T> operator+(const Vector<T>& a, const Vector<T>& b) { Vector<T> r(a); r += b; return r; }
template<class T>
inline Vector<T> operator-(const Vector<T>& a, const Vector<T>& b) { Vector<T> r(a); r -= b; return r; }

template<class T>
inline Vector<T> operator*(const Vector<T>& a, const T& s) { Vector<T> r(a); r *= s; return r; }
template<class T>
inline Vector<T> operator*(const T& s, const Vector<T>& a) { Vector<T> r(a); r *= s; return r; }
template<class T>
inline Vector<T> operator/(const Vector<T>& a, const T& s) { Vector<T> r(a); r /= s; return r; }

template<class T>
inline Vector<T> operator+(const Vector<T>& a, const T& s) { Vector<T> r(a); r += s; return r; }
template<class T>
inline Vector<T> operator-(const Vector<T>& a, const T& s) { Vector<T> r(a); r -= s; return r; }

template<class T>
inline T
product(const Vector<T> &a, const Vector<T> &b)
{
  using namespace std;

  if(a.size() != b.size())
     dmtk::err() << "** Warning: Vector sizes do not comform\n";

  typename Vector<T>::const_iterator pa = a.begin();
  typename Vector<T>::const_iterator pb = b.begin();
  int n = std::min(a.size(), b.size());
  return dot_product(n, pa, 1, pb, 1);
}

template<class T>
inline T
dot_product(const Vector<T> &a, const Vector<T> &b)
{
  return product(a,b);
}

// Streams

template<class T>
std::ostream& operator << (std::ostream& s, const Vector<T>& v)
{
  typename Vector<T>::const_iterator iter = v.begin();
  s << "(" << (*iter++);
  for(; iter != v.end(); iter++)
    s << "," << (*iter);
  s << ")";

  return s;
}


template<class T>
void
Vector<T>::write(std::ostream &s) const
{
  size_t n = this->size();
  s.write(reinterpret_cast<const char*>(&n), sizeof(size_t));

  for (T const & x : *this)
    s.write(reinterpret_cast<const char*>(&x), sizeof(T));
}

template<class T>
void
Vector<T>::read(std::istream &s)
{
  size_t n;
  s.read(reinterpret_cast<char*>(&n), sizeof(size_t));
  resize(n);

  for (T & x : *this)
    s.read(reinterpret_cast<char*>(&x), sizeof(T));
}


/////////////////////////////////////////////////////////////////

} // namespace dmtk

#endif // __DMTK_VECTOR_H__
