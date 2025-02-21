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

#ifndef __DMTK_GSLICE_ITER_H__
#define __DMTK_GSLICE_ITER_H__

// gslice iterator class. From Stroustrup's book.

#include <vector>
#include "slice_iter.h" 

namespace dmtk
{

#include "meta.h"

template<class T> class cgslice_iter;

template <class T>
class gslice_iter
{
  public:
    typedef typename std::vector<T> _V;

    gslice_iter(): v(0),s1(0,0,0),s2(0,0,0),curr1(0),curr2(0) {}
    gslice_iter(_V* vv, std::slice s1, std::slice s2):
      v(vv),s1(s1),s2(s2),curr1(0),curr2(0){};
    gslice_iter(const gslice_iter<T>& ss):
      v(ss.v),s1(ss.s1),s2(ss.s2),curr1(ss.curr1),curr2(ss.curr2){}

    gslice_iter begin()  
      {gslice_iter t = *this; t.curr1 = 0; t.curr2 = 0; return t; }
    gslice_iter end()  
      {gslice_iter t = *this; t.curr1 = t.s1.size(); t.curr2 = s2.size(); return t; }

    gslice_iter& operator=(gslice_iter<T> ss)
      { 
        size_t n1 = std::min(ss.s1.size(), s1.size());
        size_t n2 = std::min(ss.s2.size(), s2.size());
        for(size_t i = 0; i < n1; i++)
          for(size_t j = 0; j < n2; j++)
            this->operator()(i,j) = ss(i,j);
        return *this; 
      }

    gslice_iter& operator=(cgslice_iter<T> ss)
      { 
        size_t n1 = std::min(ss.size1(), s1.size());
        size_t n2 = std::min(ss.size2(), s2.size());
        for(size_t i = 0; i < n1; i++)
          for(size_t j = 0; j < n2; j++)
            this->operator()(i,j) = ss(i,j);
        return *this; 
      }

//  Operators

    template<class Expr>
    gslice_iter& operator=(const IterExpr<T,Expr>&);

    gslice_iter& operator+=(gslice_iter<T>);
    gslice_iter& operator-=(gslice_iter<T>);
    gslice_iter& operator*=(gslice_iter<T>);
    gslice_iter& operator/=(gslice_iter<T>);

    gslice_iter& operator+=(cgslice_iter<T>);
    gslice_iter& operator-=(cgslice_iter<T>);
    gslice_iter& operator*=(cgslice_iter<T>);
    gslice_iter& operator/=(cgslice_iter<T>);

    gslice_iter& operator=(_V);
    gslice_iter& operator+=(_V);
    gslice_iter& operator-=(_V);
    gslice_iter& operator*=(_V);
    gslice_iter& operator/=(_V);

    gslice_iter& operator=(const T*);
    gslice_iter& operator+=(const T*);
    gslice_iter& operator-=(const T*);
    gslice_iter& operator*=(const T*);
    gslice_iter& operator/=(const T*);

    template<class Expr>
    gslice_iter& operator+=(const IterExpr<T,Expr>&);
    template<class Expr>
    gslice_iter& operator-=(const IterExpr<T,Expr>&);
    template<class Expr>
    gslice_iter& operator*=(const IterExpr<T,Expr>&);
    template<class Expr>
    gslice_iter& operator/=(const IterExpr<T,Expr>&);

    size_t size() const { return s1.size(); }
    size_t size1() const { return s1.size(); }
    size_t size2() const { return s2.size(); }
    size_t stride1() const { return s1.stride(); }
    size_t stride2() const { return s2.stride(); }
    size_t start1() const { return s1.start(); }
    size_t start2() const { return s2.start(); }

    gslice_iter& operator++() { curr1++; return *this; }
    gslice_iter  operator++(int) {gslice_iter t = *this; curr1++; return t;}

    slice_iter<T> operator[](size_t i) { return ref(curr1 = i); }

    T& operator()(size_t i, size_t j)
      {return ref(curr1=i, curr2=j);}
    slice_iter<T>& operator*() { return ref(curr1); }

    slice_iter<T> column(size_t i) 
      { 
        return slice_iter<T>(v, std::slice(s1.start()+i*s1.stride(),s2.size(),s2.stride()));
      }
    slice_iter<T> row(size_t j) 
      { 
        return slice_iter<T>(v, std::slice(s1.start()+j*s2.start(),s1.size(),s1.stride()));
      }

    gslice_iter transpose() const
      { return gslice_iter(v, s2, s1); } 

  private:
    _V* v;
    std::slice s1;
    std::slice s2;
    size_t curr1, curr2;

    slice_iter<T> ref(size_t i) 
      { 
        return slice_iter<T>(v, std::slice(s1.start()+i*s1.stride(),s2.size(),s2.stride()));
      }

    T& ref(size_t i, size_t j) 
      { 
        return (*v)[s1.start()+i*s1.stride()+s2.start()+j*s2.stride()]; 
      }
};

template <class T>
class cgslice_iter
{
  public:
    typedef typename std::vector<T> _V;

    cgslice_iter(): v(0),s1(0,0,0),s2(0,0,0),curr1(0),curr2(0) {}
    cgslice_iter(const _V* vv, std::slice s1, std::slice s2):
      v(vv),s1(s1),s2(s2),curr1(0),curr2(0){};
    cgslice_iter(const cgslice_iter<T>& ss):
      v(ss.v),s1(ss.s1),s2(ss.s2),curr1(ss.curr1),curr2(ss.curr2){}

    cgslice_iter begin()  
      {cgslice_iter t = *this; t.curr1 = 0; t.curr2 = 0; return t; }
    cgslice_iter end()  
      {cgslice_iter t = *this; t.curr1 = s1.size(); t.curr2 = s2.size(); return t; }

/*
    cgslice_iter& operator=(const cgslice_iter<T>& ss)
      { 
        v = ss.v; s1 = ss.s1, s2 = ss.s2; curr1 = ss.curr1; curr2 = ss.curr2;
        return *this; 
      }
*/
//  Operators

    cgslice_iter& operator++() { curr1++; return *this; }
    cgslice_iter operator++(int) {cgslice_iter t = *this; curr1++; return t;}

    cslice_iter<T> operator[](size_t i) { return ref(curr1 = i); }

    T operator()(size_t i, size_t j) 
      {return ref(curr1=i, curr2=j);}

// dirty hack to get the exact address
    const T* get_pointer(size_t i, size_t j) const
      { return &(const_cast<vector<T>&>(*v)[s1.start()+i*s1.stride()+s2.start()+j*s2.stride()]); }

    cslice_iter<T>& operator*() { return ref(curr1); }

    cslice_iter<T> column(size_t i) 
      { 
        return cslice_iter<T>(v, std::slice(s1.start()+i*s1.stride(),s2.size(),s2.stride()));
      }
    cslice_iter<T> row(size_t j) 
      { 
        return cslice_iter<T>(v, std::slice(s1.start()+j*s2.start(),s1.size(),s1.stride()));
      }

    size_t size() const { return s1.size(); }
    size_t size1() const { return s1.size(); }
    size_t size2() const { return s2.size(); }
    size_t stride1() const { return s1.stride(); }
    size_t stride2() const { return s2.stride(); }
    size_t start1() const { return s1.start(); }
    size_t start2() const { return s2.start(); }

    cgslice_iter transpose() const
      { return cgslice_iter(v, s2, s1); } 

  private:
    const _V* v;
    std::slice s1;
    std::slice s2;
    size_t curr1, curr2;

    cslice_iter<T> ref(size_t i) 
      { 
        return cslice_iter<T>(v, std::slice(s1.start()+i*s1.stride(),s2.size(),s2.stride()));
      }

    T ref(size_t i, size_t j) 
      { 
        return (*v)[s1.start()+i*s1.stride()+s2.start()+j*s2.stride()]; 
      }
};


#include "gslice_implement.h"

template<class T>
template<class Expr>
inline gslice_iter<T>&
gslice_iter<T>::operator=(const IterExpr<T,Expr>& vv)
{
  for(size_t i = 0; i < s1.size(); i++)
    for(size_t j = 0; j < s2.size(); j++) 
      this->operator()(i,j) = vv(i,j); 

  return *this;
}

//////////////////////////////////////////////////////////

#define BINARY_OP(op,ap) \
template<class T> \
template<class Expr> \
inline gslice_iter<T>& \
gslice_iter<T>::op(const IterExpr<T,Expr>& mm) \
{ \
  for(size_t i = 0; i < s1.size(); i++) \
    for(size_t j = 0; j < s2.size(); j++){ \
      operator()(i,j) ap mm(i,j); \
  } \
 \
  return *this; \
}

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline gslice_iter<T>& \
gslice_iter<T>::op(gslice_iter<T> mm) \
{ \
  for(size_t i = 0; i < s1.size(); i++) \
    for(size_t j = 0; j < s2.size(); j++){ \
      operator()(i,j) ap mm(i,j); \
  } \
 \
  return *this; \
}

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline gslice_iter<T>& \
gslice_iter<T>::op(cgslice_iter<T> mm) \
{ \
  for(size_t i = 0; i < s1.size(); i++) \
    for(size_t j = 0; j < s2.size(); j++){ \
      operator()(i,j) ap mm(i,j); \
  } \
 \
  return *this; \
}

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline gslice_iter<T>& \
gslice_iter<T>::op(_V v) \
{ \
  T* _v = &v[0]; \
  for(size_t i = 0; i < s1.size(); i++) \
    for(size_t j = 0; j < s2.size(); j++){ \
      operator()(i,j) ap *_v++; \
  } \
 \
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
inline gslice_iter<T>& \
gslice_iter<T>::op(const T* _v) \
{ \
  for(size_t i = 0; i < s1.size(); i++) \
    for(size_t j = 0; j < s2.size(); j++){ \
      operator()(i,j) ap *_v++; \
  } \
 \
  return *this; \
}

BINARY_OP(operator=,=);
BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP
} // namespace dmtk 

#endif // __DMTK_GSLICE_ITER_H__


