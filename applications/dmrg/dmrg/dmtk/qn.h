#ifndef __DMTK_QN_H__
#define __DMTK_QN_H__

#include <iosfwd>
#include <alps/model.h>
#include "constants.h"
#include "bits.h"
#include "meta.h"

namespace dmtk 
{

#define QN_MAX_SIZE 4

class QN:public dmtk::Vector<alps::half_integer<short> >
{
  private:
    static Vector<std::string> _qn_name;
    static int _qn_fermion;
    static size_t QN_LAST;

    typedef alps::half_integer<short> half_integer_type;
    typedef dmtk::Vector<alps::half_integer<short> > _V;
  public:
    QN() { _V::resize(QN_MAX_SIZE); }
    QN(half_integer_type val) { _V::resize(QN_MAX_SIZE); _V::operator=(val); }
    QN(const dmtk::Vector<half_integer_type> &v) : dmtk::Vector<alps::half_integer<short> >(v) {}

    QN& operator=(const dmtk::Vector<half_integer_type>&v) { _V::operator=(v); return *this; }
    QN& operator=(const QN &_qn) { _V::operator=(_qn); return *this; }
    QN& operator=(half_integer_type  val) { _V::resize(QN_MAX_SIZE); _V::operator=(val); return *this; }

    half_integer_type operator[](size_t t) const { return _V::operator[](t); }
    half_integer_type& operator[](size_t t) { return _V::operator[](t); }
    half_integer_type operator()(size_t t) const { return _V::operator[](t); }
    half_integer_type& operator()(size_t t) { return _V::operator[](t); }

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
    static int mask(const std::string &str) { return (1 << get_qn_index(str)); }
    static void init() { _qn_name = std::string(""); _qn_fermion = 0; QN_LAST = 0; }

    // Streams

    void write(std::ostream &s) const
    {
      for(int i = 0; i < QN_LAST; i++){
        short n = _V::operator[](i).get_twice();
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
int QN::_qn_fermion = 0;

inline bool
QN::operator==(const QN &qn) const
{
  return equal(qn, QN::default_mask());
}

inline bool
QN::operator!=(const QN &qn) const
{
  QN qn1 = *this;
//  QN qn2 = qn;
  for(int i = 0; i < QN_LAST; i++){
    if(qn1[i] != qn[i]) return true;
  }
  return false;
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
    alps::half_integer<short>  i1 = qn1[i]; \
    alps::half_integer<short>  i2 = qn2[i]; \
    if(i1 == i2)  \
      continue; \
    else if(i1 ap i2)  \
      return true; \
    else; \
      return false; \
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
    alps::half_integer<short>  i1 = qn1[i]; \
    alps::half_integer<short>  i2 = qn2[i]; \
    if(i1 == i2)  \
      continue; \
    else if(i1 ap i2)  \
      return true; \
    else; \
      return false; \
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
  meta_op<QN_MAX_SIZE,DMApAdd0<alps::half_integer<short> > >::op(c, a, b);
  return c;
}

QN&
QN::operator+=(const QN &v)
{ 
  meta_op<QN_MAX_SIZE,DMApAdd0<alps::half_integer<short> > >::op(*this,*this,v); 
  return *this; 
}

QN
operator-(const QN &a, const QN &b)
{
  QN c;
  meta_op<QN_MAX_SIZE,DMApSubs0<alps::half_integer<short> > >::op(c, a, b);
  return c;
}

QN&
QN::operator-=(const QN &v)
{ 
  meta_op<QN_MAX_SIZE,DMApSubs0<alps::half_integer<short> > >::op(*this,*this,v); 
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
  return SGN(meta_and<QN_MAX_SIZE>::op(*this));  
}


} // namespace dmtk

#endif // __DMTK_QN_H__
