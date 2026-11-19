/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
#include <cstdint>
#include <vector>

namespace dmtk 
{

#define QN_MAX_SIZE 4

class QN
{
  private:
    static Vector<std::string> _qn_name;
    static int _qn_fermion;
    static int _qn_mask;

    typedef alps::half_integer<short> half_integer_type;

    alps::half_integer<short> _qn0;
    alps::half_integer<short> _qn1;
    alps::half_integer<short> _qn2;
    alps::half_integer<short> _qn3;
    alps::half_integer<short> _qn4;
    alps::half_integer<short> _qn5;

    // Writes to a quantum number the model does not define (e.g. dqn.sz() = 0
    // for a model without Sz) land here. Slots 4 and 5 are never used as
    // quantum numbers because QN_MAX_SIZE is 4, and the const accessors
    // return 0 for such names, so the stored value is never read.
    half_integer_type& unused_slot() { return _qn5; }

  public:
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
      {
        size_t idx = get_qn_index(str);
        return (idx < QN_MAX_SIZE?operator[](idx):unused_slot());
      }

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
        return (idx < QN_MAX_SIZE?operator[](idx):unused_slot()); 
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
        return (idx < QN_MAX_SIZE?operator[](idx):unused_slot()); 
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

    // The components that == compares (those in get_qn_mask()), packed into
    // one integer: equal keys mean equal quantum numbers, as long as
    // QN_LAST <= QN_MAX_SIZE. dmrg.h stores values only in the first
    // QN_MAX_SIZE components, and each is a half integer stored in 16 bits.
    uint64_t key() const
      {
        const int mask = get_qn_mask();
        uint64_t k = 0;
        for(int i = 0; i < QN_MAX_SIZE; i++)
          if(mask & (1 << i))
            k |= uint64_t(uint16_t(operator[](i).get_twice())) << (16*i);
        return k;
      }

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

// Open-addressing hash index from a key to a position in a std::vector. The
// keys are computed from the elements by a functor passed to find(). The
// index is built on first use and rebuilt when stale: after invalidate(), or
// when the vector's size or buffer, or the global QN mask, changed since it
// was built. The owner must call invalidate() when it changes an element's
// key in place. Copies keep the index; the owner passes it through copied()
// to tie it to its own vector. Not thread-safe: find() modifies the index of
// a const owner.
//
// Key needs operator== and a free function uint64_t hash_key(const Key&).
inline uint64_t hash_key(uint64_t k) { return k; }

template<class Key>
class HashIndex
{
  private:
    std::vector<Key> _keys;        // key of each element, in order
    std::vector<int> _table;       // open addressing: element position, or -1
    int _shift = 63;               // 64 - log2(_table.size())
    size_t _size = 0;              // size, buffer and QN mask it was built for
    const void *_data = 0;
    int _qn_mask = 0;
    bool _valid = false;

    static size_t slot(const Key &key, int shift)
      { return size_t((hash_key(key) * 0x9E3779B97F4A7C15ull) >> shift); }

    bool current(size_t size, const void *data) const
      { return _valid && _size == size && _data == data && _qn_mask == QN::get_qn_mask(); }

    template<class E, class KeyOf>
    void build(const std::vector<E> &v, KeyOf key_of)
      {
        const size_t n = v.size();
        int bits = 1;
        while((size_t(1) << bits) < 2*n) bits++;   // at most half full
        _keys.resize(n);
        _table.assign(size_t(1) << bits, -1);
        _shift = 64 - bits;
        const size_t mask = _table.size() - 1;
        for(size_t i = 0; i < n; i++){
          _keys[i] = key_of(v[i]);
          size_t s = slot(_keys[i], _shift);
          while(_table[s] >= 0 && !(_keys[_table[s]] == _keys[i])) s = (s + 1) & mask;
          if(_table[s] < 0) _table[s] = int(i);   // the first element of a key wins
        }
        _size = n;
        _data = v.data();
        _qn_mask = QN::get_qn_mask();
        _valid = true;
      }

  public:
    void invalidate() { _valid = false; }

    // After copying the index of vector src into the owner of vector dst:
    // keep it if it was current for src, now tied to dst.
    template<class E>
    void copied(const std::vector<E> &src, const std::vector<E> &dst)
      {
        if(current(src.size(), src.data()) && dst.size() == src.size())
          _data = dst.data();
        else
          _valid = false;
      }

    // Position of the first element e of v with key_of(e) == key, or -1.
    template<class E, class KeyOf>
    int find(const std::vector<E> &v, const Key &key, KeyOf key_of)
      {
        if(!current(v.size(), v.data())) build(v, key_of);
        const size_t mask = _table.size() - 1;
        for(size_t s = slot(key, _shift);; s = (s + 1) & mask){
          const int i = _table[s];
          if(i < 0 || _keys[i] == key) return i;
        }
      }
};

// Index from quantum number to position, for vectors whose elements have
// qn(), such as the blocks of a BMatrix or the subspaces of a PackedBasis.
class QNIndex: public HashIndex<uint64_t>
{
  private:
    struct KeyOf
    {
      template<class E>
      uint64_t operator()(const E &e) const { return e.qn().key(); }
    };

  public:
    // Position of the first element of v with quantum number qn, or -1.
    template<class E>
    int find(const std::vector<E> &v, const QN &qn)
      {
        // keys hold QN_MAX_SIZE components: compare more of them one by one
        if(QN::QN_LAST > QN_MAX_SIZE){
          for(size_t i = 0; i < v.size(); i++)
            if(v[i].qn() == qn) return int(i);
          return -1;
        }
        return HashIndex<uint64_t>::find(v, qn.key(), KeyOf());
      }
};

std::ostream& operator << (std::ostream& s, const QN& q)
{
  s << "(" << q[0];
  for(int i = 1; i < QN::QN_LAST; i++)
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
