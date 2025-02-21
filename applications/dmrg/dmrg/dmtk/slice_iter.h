/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
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

#ifndef __DMTK_SLICE_ITER_H__
#define __DMTK_SLICE_ITER_H__

// slice iterator class. From Stroustrup's book.

#include <iosfwd>
#include <vector>
#include "conj.h"
#include "meta.h"
#include "array_util.h"

namespace dmtk
{

template<class T> class cslice_iter;

template <class T>
class slice_iter
{
  public:
    typedef typename std::vector<T> _V;
    slice_iter():v(0),s(0,0,0),curr(0){};
    slice_iter(_V* vv, std::slice ss):v(vv),s(ss),curr(0){};
    slice_iter(const slice_iter<T>& ss):v(ss.v),s(ss.s),curr(ss.curr){}

    size_t size() const { return s.size(); }
    slice_iter begin() const {slice_iter t = *this; t.curr = 0; return t; }
    slice_iter end() const {slice_iter t = *this; t.curr = size(); return t; }
    size_t current() const { return curr; }

//  Asignment

    slice_iter& operator=(slice_iter<T> ss)
      { 
        int n = std:: min(ss.size(), size());
        array_copy(n, ss, *this);
//        for(int i = 0; i < n; i++) ref(i) = ss[i];
        return *this;
      }
    slice_iter& operator=(cslice_iter<T> ss)
      { 
        int n = std:: min(ss.size(), size());
        array_copy2(n, ss, *this);
//        for(int i = 0; i < n; i++) ref(i) = ss[i];
        return *this;
      }

    slice_iter& operator=(const T& v){ref(curr) = v; return *this;}

    template<class Expr>
    slice_iter& operator=(const IterExpr<T,Expr>&);

    slice_iter& operator+=(slice_iter<T>);
    slice_iter& operator-=(slice_iter<T>);
    slice_iter& operator*=(slice_iter<T>);
    slice_iter& operator/=(slice_iter<T>);

    slice_iter& operator=(_V);
    slice_iter& operator+=(_V);
    slice_iter& operator-=(_V);
    slice_iter& operator*=(_V);
    slice_iter& operator/=(_V);

    slice_iter& operator=(const T*);
    slice_iter& operator+=(const T*);
    slice_iter& operator-=(const T*);
    slice_iter& operator*=(const T*);
    slice_iter& operator/=(const T*);

    template<class Expr>
    slice_iter& operator+=(const IterExpr<T,Expr>&);
    template<class Expr>
    slice_iter& operator-=(const IterExpr<T,Expr>&);
    template<class Expr>
    slice_iter& operator*=(const IterExpr<T,Expr>&);
    template<class Expr>
    slice_iter& operator/=(const IterExpr<T,Expr>&);

//  Operators

    slice_iter& operator++() { curr++; return *this; }
    slice_iter  operator++(int) {slice_iter t = *this; curr++; return t;}
    slice_iter& operator--() { curr--; return *this; }
    slice_iter  operator--(int) {slice_iter t = *this; curr--; return t;}
    T& operator[](size_t i){return ref(curr=i);}
    const T& operator[](size_t i) const {return ref(i);}
    T& operator()(size_t i){return ref(curr=i);}
    const T& operator()(size_t i) const {return ref(i);}
    T& operator*() { return ref(curr);}
    T& operator*() const { return ref(curr);}

    bool operator==(const slice_iter<T>& q)
      {return curr==q.curr && s.stride()==q.s.stride() && s.start() == q.s.start();}

    bool operator!=(const slice_iter<T>& q){return !(*this==q);}

    bool operator<(const slice_iter<T>& q)
      {return curr<q.curr && s.stride()==q.s.stride() && s.start() == q.s.start();}

    bool operator>(const slice_iter<T>& q)
      {return curr>q.curr && s.stride()==q.s.stride() && s.start() == q.s.start();}

//  IterExpr auxiliary methods

    size_t size1() const { return 1; }
    size_t size2() const { return s.size(); }

  private:
    _V* v;
    std::slice s;
    size_t curr;
    T& ref(size_t i) const { return (*v)[s.start()+i*s.stride()];}
};

template <class T>
class cslice_iter
{
  public:
    typedef typename std::vector<T> _V;

    cslice_iter():v(0),s(0,0,0),curr(0){};
    cslice_iter(const _V* vv, std::slice ss):v(vv),s(ss),curr(0){};
    cslice_iter(const cslice_iter<T>& ss):v(ss.v),s(ss.s),curr(ss.curr){}

    size_t size() const { return s.size(); }
    cslice_iter begin() const {cslice_iter t = *this; t.curr = 0; return t; }
    cslice_iter end() const {cslice_iter t = *this; t.curr = size(); return t; }
    size_t current() const { return curr; }

//  Operators

    const cslice_iter& operator++() { curr++; return *this; }
    const cslice_iter operator++(int) {cslice_iter t = *this; curr++; return t;}
    T operator[](size_t i) {return ref(curr=i);}
    T operator()(size_t i) {return ref(curr=i);}
    T operator*() { return ref(curr);}
    T operator*() const { return ref(curr);}

// Dirty hack to get the exact address
    const T* get_pointer(size_t i) const
      { return &(const_cast<vector<T>&>(*v)[s.start()+i*s.stride()]); }


    bool operator==(const cslice_iter<T>& q)
      {return curr==q.curr && s.stride()==q.s.stride() && s.start() == q.s.start();}

    bool operator!=(const cslice_iter<T>& q){return !(*this==q);}

    bool operator<(const cslice_iter<T>& q)
      {return curr < q.curr && s.stride()==q.s.stride() && s.start() == q.s.start();}

    bool operator>(const cslice_iter<T>& q)
      {return curr > q.curr && s.stride()==q.s.stride() && s.start() == q.s.start();}

//  IterExpr auxiliary methods

    size_t size1() const { return 1; }
    size_t size2() const { return s.size(); }

  private:
    const _V* v;
    std::slice s;
    size_t curr;
    T ref(size_t i) const { return (*v)[s.start()+i*s.stride()];}
};


#include "slice_implement.h"

template<class T>
template<class Expr>
inline slice_iter<T>&
slice_iter<T>::operator=(const IterExpr<T,Expr>& vv)
{
  curr = 0;
  array_copy2(vv.size2(), vv, *this);
//  for(uint i = 0; i < size(); i++) { *this->operator++(0) = vv[i]; }
  curr = 0;
  return *this;
}

#define BINARY_OP(op,ap) \
template<class T> \
template<class Expr> \
inline slice_iter<T>& \
slice_iter<T>::op(const IterExpr<T,Expr>& vv) \
{ \
  curr = 0; \
  for(int i = 0; i < size(); i++) { *this->operator++(0) ap vv[i]; } \
  curr = 0; \
  return *this; \
} 

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline slice_iter<T>& \
slice_iter<T>::op(slice_iter<T> vv) \
{ \
  curr = 0; \
  for(int i = 0; i < size(); i++) { *this->operator++(0) ap vv[i]; } \
  curr = 0; \
  return *this; \
} 

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline slice_iter<T>& \
slice_iter<T>::op(_V vv) \
{ \
  curr = 0; \
  int n = size(); \
  T* _v = &vv[0]; \
  while(n--) { *this->operator++(0) ap *_v++; } \
  curr = 0; \
  return *this; \
} 

BINARY_OP(operator=,=);
BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline slice_iter<T>& \
slice_iter<T>::op(const T* _v) \
{ \
  curr = 0; \
  int n = size(); \
  while(n--) { *this->operator++(0) ap *_v++; } \
  curr = 0; \
  return *this; \
} 

BINARY_OP(operator=,=);
BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

template<class T>
inline T
product(cslice_iter<T> a, cslice_iter<T> b)
{
  if(a.size() != b.size())
     cerr << "** Warning: sizes do not comform\n";

  T r = T(0);

  for(size_t i = 0; i < std::min(a.size(), b.size()); i++) 
    r += std::conj(a[i])*b[i]; 

  return r;
}

template<class T>
inline T
product(slice_iter<T> a, slice_iter<T> b)
{
  if(a.size() != b.size())
     cerr << "** Warning: sizes do not comform\n";

  T r = T(0);

  for(size_t i = 0; i < std::min(a.size(), b.size()); i++) 
    r += std::conj(a[i])*b[i]; 

  return r;
}

} // namespace dmtk 

#endif // __DMTK_SLICE_ITER_H__


