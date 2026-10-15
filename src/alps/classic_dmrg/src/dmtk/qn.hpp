// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/qn.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#ifndef __DMTK_QN_H__
#define __DMTK_QN_H__

#include <iosfwd>
#include "dmtk/half_integer.hpp"
#include "dmtk/constants.hpp"
#include "dmtk/bits.hpp"
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

    using half_integer_type = dmtk::half_integer<short>;

    dmtk::half_integer<short> _qn0;
    dmtk::half_integer<short> _qn1;
    dmtk::half_integer<short> _qn2;
    dmtk::half_integer<short> _qn3;
    dmtk::half_integer<short> _qn4;
    dmtk::half_integer<short> _qn5;

  public:
    // Private engine metadata must outlive every QN-dependent object. The
    // invocation swaps the complete registry before construction and after
    // destruction, including the mask intentionally untouched by legacy init().
    struct Registry {
      Vector<std::string> names = Vector<std::string>(4);
      int fermion = 0;
      int mask = 0;
      size_t count = 0;
    };

    static void swap_registry(Registry& other) noexcept
    {
      _qn_name.swap(other.names);
      std::swap(_qn_fermion, other.fermion);
      std::swap(_qn_mask, other.mask);
      std::swap(QN_LAST, other.count);
    }

    static size_t QN_LAST;
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
        s.write(reinterpret_cast<const char*>(&n), sizeof(short));
      }
    }

    void read(std::istream &s)
    {
      for(int i = 0; i < QN_LAST; i++){
        short n;
        s.read(reinterpret_cast<char*>(&n), sizeof(short));
        operator[](i).set_half(n);
      }
    }

};

// C++17 inline member-variable definitions. Without `inline`, every TU
// that includes qn.hpp gets its own copy of these static members — fine
// for the original single-TU dmrg binary, broken when dmrg_impl is
// also linked into a alps nanobind extension. Multiple definitions
// in the same address space caused duplicated startup init and a
// double-destroy SIGABRT at process exit.
inline size_t QN::QN_LAST = 0;
inline dmtk::Vector<std::string> QN::_qn_name = dmtk::Vector<std::string>(4);
inline int QN::_qn_mask = 0;
inline int QN::_qn_fermion = 0;

inline std::ostream& operator << (std::ostream& s, const QN& q)
{
  s << "(" << q[0];
  for(size_t i = 1; i < QN::QN_LAST; i++)
    s << "," << q[i];
  s << ")";

  return s;
}

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
      dmtk::half_integer<short>  i1 = qn1[i]; \
      dmtk::half_integer<short>  i2 = qn2[i]; \
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
      dmtk::half_integer<short>  i1 = qn1[i]; \
      dmtk::half_integer<short>  i2 = qn2[i]; \
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

inline QN
operator+(const QN &a, const QN &b)
{
  QN c;
  for (int i = 0; i < QN_MAX_SIZE; ++i)
    c[i] = a[i] + b[i];
  return c;
}

inline QN&
QN::operator+=(const QN &v)
{
  for (int i = 0; i < QN_MAX_SIZE; ++i)
    (*this)[i] += v[i];
  return *this;
}

inline QN
operator-(const QN &a, const QN &b)
{
  QN c;
  for (int i = 0; i < QN_MAX_SIZE; ++i)
    c[i] = a[i] - b[i];
  return c;
}

inline QN&
QN::operator-=(const QN &v)
{
  for (int i = 0; i < QN_MAX_SIZE; ++i)
    (*this)[i] -= v[i];
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
