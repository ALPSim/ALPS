/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003 by Matthias Troyer <troyer@comp-phys.org>
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

#ifndef ALPS_APPLICATIONS_MC_SPIN_TINYVEC_H_
#define ALPS_APPLICATIONS_MC_SPIN_TINYVEC_H_

#include <algorithm>
#include <iostream>
#include "helper.h"
#include "matrices.h"

template <class T, int N>
struct ScaledTinyVector;

template <class T, int N>
class TinyVector {
public:
  TinyVector(T x=0.) { std::fill(vec_,vec_+N,x); vec_[0]=x; }
  TinyVector(T x, T y) { std::fill(vec_,vec_+N,x); vec_[0]=x; vec_[1]=y; }
  TinyVector(T x, T y, T z) { std::fill(vec_,vec_+N,x); vec_[0]=x;vec_[1]=y; vec_[2]=z; }
  TinyVector(std::string&);
  TinyVector(std::string&, alps::Parameters& parms);
  inline void operator*=(T y);
  void operator/=(T y) { operator*=(1./y);}
  template <class Y> inline void operator-=(Y y);
  template <class Y> inline void operator+=(Y y);
  T operator[](int i) const { return vec_[i];}
  T& operator[](int i) { return vec_[i];}
  T get_length2();

  TinyVector<T,N> mult(const SquareMatrix<T,N> m) {
    TinyVector<T,N> res;
    for (int i=0; i<N; i++) {
      res[i] = 0;
      for (int j=0; j<N;j++)
        res[i] += vec_[i]*m.getElement_nc(i,j);
    }
    return res;
  }

  TinyVector<T,N> operator* (const SquareMatrix<T,N> m) {
    TinyVector<T,N> res;
    for (int i=0; i<N; i++) {
      res[i] = 0;
      for (int j=0; j<N; j++)
        res[i]+= vec_[i]*m.getElement_nc(i,j);
    }
    return res;
  }
 
  TinyVector<T,N> operator* (const T m) {
    TinyVector<T,N> res;
    for (int i=0; i<N; i++)
      res[i] = vec_[i]*m;
    return res;
  } 

  TinyVector<T,N> operator-(const TinyVector<T,N> v2) const {
    TinyVector<T,N> res;
    for (int i=0; i<N; i++)
      res[i] = vec_[i] - v2[i];
    return res;
  }
  friend alps::ODump& operator << (alps::ODump& dump, const TinyVector<T,N>& m) 
  { for (int i=0;i<N;++i) dump << m[i]; return dump;}
  friend std::ostream& operator << (std::ostream& out, const TinyVector<T,N>& m)
  { for (int i=0;i<N;++i) out << m[i] << " "; return out;}
  friend alps::IDump& operator >> (alps::IDump& dump, TinyVector<T,N>& m) 
  { for (int i=0;i<N;++i) dump >> m[i]; return dump;}
private:
  T vec_[N];
};

/**
 * Reads in a string and computes thereof the elements for the vector.
 */
template<class T, int N>
TinyVector<T,N>::TinyVector(std::string& str) 
{
  Helper<T> myHelper(str);
  int str_count = myHelper.elemcount();
  if (str_count == 1) {
    T x = myHelper.getNextElement();
    std::fill(vec_,vec_+N,0.);
    vec_[N-1] = x;
  } else if (str_count == N)
    for (int j=0; j<N; j++)
      vec_[j] = myHelper.getNextElement();
  else {
    std::cerr << "too few input values for the diagonal matrix: \n1 or"
              << N << " values are needed, found " << str_count << ".\n"
              << "input string was " << str << ".\n";
    boost::throw_exception(std::runtime_error("too few input values"));
  }
  std::cerr << "vector set to " << *this << "\n";
}  

template<class T, int N>
TinyVector<T,N>::TinyVector(std::string& str, alps::Parameters& parms)
{
  Helper<T> myHelper(str,parms);
  int str_count = myHelper.elemcount();
  if (str_count == 1) {
    T x = myHelper.getNextElement();
    std::fill(vec_,vec_+N,0.);
    vec_[N-1] = x;
  } else if (str_count == N)
    for (int j=0; j<N; j++)
      vec_[j] = myHelper.getNextElement();
  else {
    std::cerr << "too few input values for the diagonal matrix: \n1 or"
              << N << " values are needed, found " << str_count << ".\n"
              << "input string was " << str << ".\n";
    boost::throw_exception(std::runtime_error("too few input values"));
  }
}

template <class T, int N, int I>
struct meta {
  static T dot(const TinyVector<T,N>& x, const TinyVector<T,N>& y) {
    return meta<T,N,I-1>::dot(x,y) + x[I-1]*y[I-1];}
  static void scale(TinyVector<T,N>& x, T y) {
    x[I-1]*=y; meta<T,N,I-1>::scale(x,y);
   }
  template <class Y>
  static void subtract(TinyVector<T,N>& x, Y y) {
    x[I-1]-=y[I-1];
    meta<T,N,I-1>::subtract(x,y);
  }
  template <class Y>
  static void add(TinyVector<T,N>& x, Y y) {
    x[I-1]+=y[I-1];
    meta<T,N,I-1>::add(x,y);
  }
  static T length2(const TinyVector<T,N>& x) {
    return meta<T,N,I-1>::length2(x) + x[I-1]*x[I-1]; }
};

template <class T, int N>
struct meta<T,N,0> {
  static T dot(const TinyVector<T,N>& , const TinyVector<T,N>& ) {return 0;}
  static void scale(TinyVector<T,N>& , T ) {}
  template <class Y>
  static void subtract(TinyVector<T,N>& , Y ) {}
  template <class Y>
  static void add(TinyVector<T,N>& , Y ) {}
  static T length2(const TinyVector<T,N>&) { return 0.; }
};

template <class T, int N>
struct ScaledTinyVector {
  ScaledTinyVector(T x, const TinyVector<T,N>& y) : scale(x),vec(y) {}
  T scale;
  const TinyVector<T,N>& vec;
  T operator[](int i) { return scale*vec[i];}
};

template <class T, int N> template<class Y>
inline void TinyVector<T,N>::operator-=(Y y)
{
  meta<T,N,N>::subtract(*this,y);
}

template <class T, int N> template<class Y>
inline void TinyVector<T,N>::operator+=(Y y)
{
  meta<T,N,N>::add(*this,y);
}

template <class T, int N>
inline T dot (const TinyVector<T,N>& x, const TinyVector<T,N>& y)
{
  return meta<T,N,N>::dot(x,y);
}

template <class T, int N>
inline void TinyVector<T,N>::operator*=(T x)
{
  meta<T,N,N>::scale(*this,x);
}

template <class T, int N>
inline T TinyVector<T,N>::get_length2()
{ return meta<T,N,N>::length2(*this); }

template <class T, int N>
inline ScaledTinyVector<T,N> operator* (T x, const TinyVector<T,N>& y)
{
  return ScaledTinyVector<T,N>(x,y);
}

template <class T, int N>
inline ScaledTinyVector<T,N> operator* (const TinyVector<T,N>& y, T x)
{
  return ScaledTinyVector<T,N>(x,y);
}

namespace std {
template <class T, int N>
inline T abs (const TinyVector<T,N>& x)
{
  return std::sqrt(dot(x,x));
}
}

template <class T, int N>
std::ostream & operator<<  (std::ostream& os, const TinyVector<T,N>& x)
{
  os << "(";
  for (int i=0;i<N;++i)
    os << x[i] << " ";
  os << ")";
  return os;
}

#endif
