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

#ifndef __DMTK_QN_H__
#define __DMTK_QN_H__

#include <iosfwd>
#include <alps/model.h>
#include "constants.h"
#include "bits.h"
#include "meta.h"
#include <cassert>

namespace dmtk 
{

#define QN_MAX_SIZE 4

class QN
{
  private:
    static Vector<std::string> _qn_name;
    static int _qn_fermion;
    static int _qn_mask;
    static size_t QN_LAST;

    typedef alps::half_integer<short> half_integer_type;

    alps::half_integer<short> _qn0;
    alps::half_integer<short> _qn1;
    alps::half_integer<short> _qn2;
    alps::half_integer<short> _qn3;
    alps::half_integer<short> _qn4;
    alps::half_integer<short> _qn5;

  public:
    QN(): _qn0(0),_qn1(0),_qn2(0),_qn3(0),_qn4(0),_qn5(0) {}
    QN(half_integer_type val): _qn0(val),_qn1(val),_qn2(val),_qn3(val),_qn4(val),_qn5(val) {}
    QN(const QN &_qn):
      _qn0(_qn._qn0),_qn1(_qn._qn1),_qn2(_qn._qn2),_qn3(_qn._qn3),_qn4(_qn._qn4),_qn5(_qn._qn5) {}


    QN& operator=(const QN &_qn)
      { _qn0 = _qn._qn0; _qn1 = _qn._qn1; _qn2 = _qn._qn2; _qn3 = _qn._qn3; _qn4 = _qn._qn4; _qn5 = _qn._qn5; return *this; }
    QN& operator=(half_integer_type  val) { _qn0=_qn1=_qn2=_qn3=_qn4=_qn5=val; return *this; }

    half_integer_type operator[](int i) const {
      switch(i){
        case 0:
          return _qn0;
        case 1:
          return _qn1;
        case 2:
          return _qn2;
        case 3:
          return _qn3;
        case 4:
          return _qn4;
        case 5:
          return _qn5;
        default:
          return DMTK_ERROR;
       }
    }
    half_integer_type & operator[](int i) {
      assert(i >= 0);
      assert(i <= 5);
      switch(i){
        case 0:
          return _qn0;
        case 1:
          return _qn1;
        case 2:
          return _qn2;
        case 3:
          return _qn3;
        case 4:
          return _qn4;
        case 5:
        default:
          return _qn5;
       }
    }



    half_integer_type operator[](const std::string& str) const
      { return operator[](get_qn_index(str)); }
    half_integer_type& operator[](const std::string& str)
      { return operator[](get_qn_index(str)); }

    bool equal(const QN &qn, int mask) const;
    bool operator==(const QN &qn) const;
    bool operator!=(const QN &qn) const;
    bool operator>(const QN &qn) const;
    bool operator<(const QN &qn) const;
    bool operator>=(const QN &qn) const;
    bool operator<=(const QN &qn) const;

    QN& operator+=(const QN &v);
    QN& operator-=(const QN &v);

    half_integer_type sz() const 
      { 
        size_t idx = get_qn_index("Sz"); 
        if(idx >= QN_MAX_SIZE) idx = get_qn_index("SZ"); 
        if(idx >= QN_MAX_SIZE) idx = get_qn_index("sz"); 
        return (idx < QN_MAX_SIZE?operator[](idx):0); 
      }
    half_integer_type& sz() 
      { 
        size_t idx = get_qn_index("Sz"); 
        if(idx >= QN_MAX_SIZE) idx = get_qn_index("SZ"); 
        if(idx >= QN_MAX_SIZE) idx = get_qn_index("sz"); 
        return (operator[](idx)); 
      }
    half_integer_type n() const 
      { 
        size_t idx = get_qn_index("N"); 
        if(idx >= QN_MAX_SIZE) idx = get_qn_index("n"); 
        return (idx < QN_MAX_SIZE?operator[](idx):0); 
      }
    half_integer_type& n() 
      { 
        size_t idx = get_qn_index("N"); 
        if(idx >= QN_MAX_SIZE) idx = get_qn_index("n"); 
        return (operator[](idx)); 
      }

    int fermion_sign() const;
    static size_t max_index() { return QN::QN_LAST; }
    static const std::string& qn_name(size_t i) { return _qn_name[i]; }
        
 
    static size_t get_qn_index(const std::string& str) 
      {
        Vector<std::string>::iterator iter;
        size_t idx = 0;
        for(iter = _qn_name.begin(); iter != _qn_name.end(); iter++, idx++){
          if(*iter == str) return idx;
        }
        return 999;
      };

    static void set_qn_index(size_t idx, const std::string& str)
      { _qn_name(idx) = str; }

    static size_t add_qn_index(const std::string& str, bool _fermion = false)
      {
         size_t idx = get_qn_index(str);
         if(idx == 999) {
           idx = QN_LAST;
           _qn_name.resize(QN_LAST+1);
           set_qn_index(QN_LAST, str);
           if(_fermion) _qn_fermion |= (1 << QN_LAST);
           QN_LAST++;
         }
         return idx;
      } 

    static int qn_index_fermion(size_t idx) { return (IBITS(_qn_fermion,idx)); } 

    static int default_mask() { return ((1 << QN_LAST)-1); }
    static int get_qn_mask() { return (_qn_mask); }
    static int mask(const std::string &str) { return (1 << get_qn_index(str)); }
    static void init() { _qn_name = std::string(""); _qn_fermion = 0; QN_LAST = 0; }

    static void set_qn_mask(int qn_mask)
    {
     _qn_mask = qn_mask;
    }

    // Streams

    void write(std::ostream &s) const
    {
      for(int i = 0; i < QN_LAST; i++){
        short n = this->operator[](i).get_twice();
        s.write((const char *)&n, sizeof(short));
      }
    }

    void read(std::istream &s)
    {
      for(int i = 0; i < QN_LAST; i++){
        short n;
        s.read((char *)&n, sizeof(short));
        operator[](i).set_half(n);
      }
    }

};

size_t QN::QN_LAST = 0;
dmtk::Vector<std::string> QN::_qn_name = dmtk::Vector<std::string>(4);
int QN::_qn_mask = 0;
int QN::_qn_fermion = 0;

inline bool
QN::operator==(const QN &qn) const
{
  return equal(qn, QN::get_qn_mask());
}

inline bool
QN::operator!=(const QN &qn) const
{
  return (!equal(qn, QN::get_qn_mask()));
}

inline bool
QN::equal(const QN &qn, int mask) const
{
  QN qn1 = *this;
//  QN qn2 = qn;
  for(int i = 0; i < QN_LAST; i++){
    if((mask & (1 << i)) && qn1[i] != qn[i]) return false;
  }
  return true;
}

#define OP_EXCLUSIVE(op,ap) \
inline bool \
op(const QN &qn) const \
{ \
  QN qn1 = *this; \
  QN qn2 = qn; \
  for(int i = 0; i < QN_LAST; i++){ \
    if(_qn_mask & (1 << i)){ \
      alps::half_integer<short>  i1 = qn1[i]; \
      alps::half_integer<short>  i2 = qn2[i]; \
      if(i1 == i2)  \
        continue; \
      else if(i1 ap i2)  \
        return true; \
      else; \
        return false; \
    }\
  } \
  \
  return false; \
} 

OP_EXCLUSIVE(QN::operator>,>)
OP_EXCLUSIVE(QN::operator<,<)
#undef OP_EXCLUSIVE

#define OP_INCLUSIVE(op,ap) \
inline bool \
op(const QN &qn) const \
{ \
  QN qn1 = *this; \
  QN qn2 = qn; \
  for(int i = 0; i < QN_LAST; i++){ \
    if(_qn_mask & (1 << i)){ \
      alps::half_integer<short>  i1 = qn1[i]; \
      alps::half_integer<short>  i2 = qn2[i]; \
      if(i1 == i2)  \
        continue; \
      else if(i1 ap i2)  \
        return true; \
      else; \
        return false; \
    }\
  } \
  \
  return true; \
} 

OP_INCLUSIVE(QN::operator>=,>=)
OP_INCLUSIVE(QN::operator<=,<=)
#undef OP_INCLUSIVE

template <int I, class BinOp> 
struct meta_op{
  static void op(QN& a, const QN& b, const QN &c)
    { a[I] = BinOp::apply(b[I],c[I]); meta_op<I-1,BinOp>::op(a,b,c); }
};

template<class BinOp>
struct meta_op<0,BinOp>{
  static void op(QN& a, const QN& b, const QN &c)
    { a[0] = BinOp::apply(b[0],c[0]); }
};

QN
operator+(const QN &a, const QN &b)
{
  QN c;
  meta_op<QN_MAX_SIZE-1,DMApAdd0<alps::half_integer<short> > >::op(c, a, b);
  return c;
}

QN&
QN::operator+=(const QN &v)
{ 
  meta_op<QN_MAX_SIZE-1,DMApAdd0<alps::half_integer<short> > >::op(*this,*this,v); 
  return *this; 
}

QN
operator-(const QN &a, const QN &b)
{
  QN c;
  meta_op<QN_MAX_SIZE-1,DMApSubs0<alps::half_integer<short> > >::op(c, a, b);
  return c;
}

QN&
QN::operator-=(const QN &v)
{ 
  meta_op<QN_MAX_SIZE-1,DMApSubs0<alps::half_integer<short> > >::op(*this,*this,v); 
  return *this; 
}

template <int I> 
struct meta_and{
  static int op(const QN& a)
    { return (QN::qn_index_fermion(I)*a[I].get_twice()/2 + meta_and<I-1>::op(a)); }
};

template<>
struct meta_and<0>{
  static int op(const QN& a)
    { return QN::qn_index_fermion(0)*a[0].get_twice()/2; }
};

inline int 
QN::fermion_sign() const
{
  return SGN(meta_and<QN_MAX_SIZE-1>::op(*this));  
}


} // namespace dmtk

#endif // __DMTK_QN_H__
