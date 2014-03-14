/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef __DMTK_STATE_H__
#define __DMTK_STATE_H__

#include <iosfwd>
#include <vector>
#include "enums.h"
#include "range.h"
#include "slice_iter.h"
#include "gslice_iter.h"
#include "qn.h"
#include "subspace.h"
#include "basis.h"
#include "state_slice.h"
#include "block_matrix.h"
#include "util.h"
#include "lanczos.cc"

using namespace std;

namespace dmtk
{

#include "meta.h"

template <class T> class VectorState;
template<class T>
void new_seed(const VectorState<T>& v, VectorState<T>& res,
              const BMatrix<T>& rho1, const BMatrix<T> &rho2,
              const Basis& basis1, const Basis& basis2, int position);

class StateSpace
{
  private:
    SubSpace _s1;
    SubSpace _s2;
    SubSpace _s3;
    SubSpace _s4;
    int _start;  /* negative means error */

    inline const SubSpace& ref(size_t i) const
      { 
        switch(i){
          case 1: return _s1;
          case 2: return _s2;
          case 3: return _s3;
          case 4: return _s4;
          default: return _s1;
        }
      }
    inline SubSpace& ref(size_t i)
      { 
        switch(i){
          case 1: return _s1;
          case 2: return _s2;
          case 3: return _s3;
          case 4: return _s4;
          default: return _s1;
        }
      }

  public:

    StateSpace(): _start(0) {}
    StateSpace(SubSpace s1, SubSpace s2, SubSpace s3, SubSpace s4, size_t start = 0) : _start(start) { _s1 = s1; _s2 = s2; _s3 = s3; _s4 = s4; } 

    StateSpace& operator=(const StateSpace &s)
      { _start = s._start; _s1 = s._s1; _s2 = s._s2; _s3 = s._s3; _s4 = s._s4; return *this; }

    SubSpace& operator()(size_t i) { return ref(i); } 
    SubSpace operator()(size_t i) const { return ref(i); }
    SubSpace& operator[](size_t i) { return ref(i); } 
    SubSpace operator[](size_t i) const { return ref(i); }

    bool operator==(const StateSpace&) const;
    bool operator!=(const StateSpace&) const;
    bool operator>(const StateSpace&) const;
    bool operator>=(const StateSpace&) const;
    bool operator<(const StateSpace&) const;
    bool operator<=(const StateSpace&) const;

    size_t dim() const { return _s1.dim() * _s2.dim() * _s3.dim() * _s4.dim(); }
    int start() const { return _start; } /* negative means error */

    template<class T>
    friend class VectorState;
}; 

inline bool
StateSpace::operator==(const StateSpace &s) const
{
  for(int i = 1; i < 5; i++)
    if(ref(i).qn() != s.ref(i).qn()) return false;
  return true;
}

inline bool
StateSpace::operator!=(const StateSpace &s) const
{
  for(int i = 1; i < 5; i++)
    if(ref(i).qn() != s.ref(i).qn()) return true;
  return false;
}

#define OP_EXCLUSIVE(op,ap) \
inline bool \
op(const StateSpace &s) const \
{ \
  for(int i = 1; i < 5; i++){ \
    QN i1 = ref(i).qn(); \
    QN i2 = s.ref(i).qn(); \
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

OP_EXCLUSIVE(StateSpace::operator>,>)
OP_EXCLUSIVE(StateSpace::operator<,<)
#undef OP_EXCLUSIVE

#define OP_INCLUSIVE(op,ap) \
inline bool \
op(const StateSpace &s) const \
{ \
  for(int i = 1; i < 5; i++){ \
    QN i1 = ref(i).qn(); \
    QN i2 = s.ref(i).qn(); \
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

OP_INCLUSIVE(StateSpace::operator>=,>=)
OP_INCLUSIVE(StateSpace::operator<=,<=)
#undef OP_INCLUSIVE

template <class T>
class VectorState: public Vector<T>
{
  private:
    PackedBasis _b1;
    PackedBasis _b2;
    PackedBasis _b3;
    PackedBasis _b4;
    QN _qn;
    int qn_constrained;
    std::vector<StateSpace> qn_space;
    Vector<size_t> index; // for condensed vectors

    void resize_constrained();
    void copy_constrained(const Vector<T>& v);
    T& ref_constrained(size_t i1, size_t i2, size_t i3, size_t i4);
    T ref_constrained(size_t i1, size_t i2, size_t i3, size_t i4) const;

    VectorState condense1(int mask) const;
    VectorState decondense1(int mask, const VectorState<T> &orig) const;

    void init(const Basis &b1, const Basis &b2, 
              const Basis &b3, const Basis &b4)
      { _b1 = b1.subspaces(); _b2 = b2.subspaces(); _b3 = b3.subspaces(); _b4 = b4.subspaces(); }
    void init(const PackedBasis &b1, const PackedBasis &b2, 
              const PackedBasis &b3, const PackedBasis &b4)
      { _b1 = b1; _b2 = b2; _b3 = b3; _b4 = b4; }

  public:

    typedef std::vector<StateSpace> _V;
    typedef std::vector<StateSpace>::iterator iterator;
    typedef std::vector<StateSpace>::const_iterator const_iterator;

    VectorState(): qn_constrained(0) {}
    VectorState(size_t _d): Vector<T>(_d), qn_constrained(0){}
    VectorState(size_t _d1, size_t _d2, size_t _d3, size_t _d4): Vector<T>(_d1*_d2*_d3*_d4), qn_constrained(0)
        { init(Basis(_d1).pack(),Basis(_d2).pack(),Basis(_d3).pack(),Basis(_d4).pack()); }
    VectorState(const Basis &b1, const Basis &b2, 
                const Basis &b3, const Basis &b4):
      Vector<T>(b1.dim()*b2.dim()*b3.dim()*b4.dim()), qn_constrained(0) 
        { init(b1,b2,b3,b4); }
    VectorState(const Vector<T> &v, 
                const Basis &b1, const Basis &b2, 
                const Basis &b3, const Basis &b4):
      Vector<T>(v), qn_constrained(0) { init(b1,b2,b3,b4); }
    VectorState(const Basis &b1, const Basis &b2, 
                const Basis &b3, const Basis &b4, const QN &qn):
      _qn(qn), qn_constrained(0)
      { init(b1,b2,b3,b4); resize_constrained(); }
    VectorState(const PackedBasis &b1, const PackedBasis &b2, 
                const PackedBasis &b3, const PackedBasis &b4):
      Vector<T>(b1.dim()*b2.dim()*b3.dim()*b4.dim()), qn_constrained(0) 
        { init(b1,b2,b3,b4); }
    VectorState(const Vector<T> &v, 
                const PackedBasis &b1, const PackedBasis &b2, 
                const PackedBasis &b3, const PackedBasis &b4):
      Vector<T>(v), qn_constrained(0) { init(b1,b2,b3,b4); }
    VectorState(const PackedBasis &b1, const PackedBasis &b2, 
                const PackedBasis &b3, const PackedBasis &b4, const QN &qn,
                int qn_mask = 0):
      _qn(qn), qn_constrained(qn_mask)
      { init(b1,b2,b3,b4); resize_constrained(); }

    VectorState(const Vector<T> &v, 
                const Basis &b1, const Basis &b2, 
                const Basis &b3, const Basis &b4, const QN &qn,
                int qn_mask = 0):
      _qn(qn), qn_constrained(qn_mask)
      { init(b1,b2,b3,b4); resize_constrained(); copy_constrained(v); }
    VectorState(const Vector<T> &v, 
                const PackedBasis &b1, const PackedBasis &b2, 
                const PackedBasis &b3, const PackedBasis &b4, const QN &qn,
                int qn_mask = 0):
      _qn(qn), qn_constrained(qn_mask)
      { init(b1,b2,b3,b4); resize_constrained(); copy_constrained(v); }

    VectorState(const VectorState<T>& v):
      Vector<T>(v), _b1(v._b1), _b2(v._b2), _b3(v._b3), _b4(v._b4), _qn(v._qn), qn_constrained(v.qn_constrained), qn_space(v.qn_space), index(v.index) {}


    VectorState& operator=(const T& v)
      { Vector<T>::operator=(v); return *this; }
    
    VectorState& operator=(const VectorState<T>& v)
      { Vector<T>::operator=(v); _b1 = v._b1; _b2 = v._b2; _b3 = v._b3; _b4 = v._b4; _qn = v._qn; qn_constrained = v.qn_constrained; qn_space = v.qn_space; index = v.index; return *this; }

    VectorState& operator=(const Vector<T>& v)
      {
        if(qn_constrained == 0)
          Vector<T>::operator=(v);
        else 
          copy_constrained(v);

        return *this;
      };

    template<class Expr>
    VectorState& operator=(const IterExpr<T,Expr>& exp)
      { Vector<T>::operator=(exp); return *this; }

    VectorState& operator+=(const VectorState<T> &v)
      { Vector<T>::operator+=(v); return *this; }
    VectorState& operator-=(const VectorState<T> &v)
      { Vector<T>::operator-=(v); return *this; }
    VectorState& operator*=(const VectorState<T> &v)
      { Vector<T>::operator*=(v); return *this; }
    VectorState& operator/=(const VectorState<T> &v)
      { Vector<T>::operator/=(v); return *this; }

    VectorState& operator+=(const T &v)
      { Vector<T>::operator+=(v); return *this; }
    VectorState& operator-=(const T &v)
      { Vector<T>::operator-=(v); return *this; }
    VectorState& operator*=(const T &v)
      { Vector<T>::operator*=(v); return *this; }
    VectorState& operator/=(const T &v)
      { Vector<T>::operator/=(v); return *this; }

    template<class Expr>
    VectorState& operator+=(const IterExpr<T,Expr>&exp)
      { Vector<T>::operator+=(exp); return *this; }
    template<class Expr>
    VectorState& operator-=(const IterExpr<T,Expr>&exp)
      { Vector<T>::operator-=(exp); return *this; }
    template<class Expr>
    VectorState& operator*=(const IterExpr<T,Expr>&exp)
      { Vector<T>::operator*=(exp); return *this; }
    template<class Expr>
    VectorState& operator/=(const IterExpr<T,Expr>&exp)
      { Vector<T>::operator/=(exp); return *this; }

    T& operator[](size_t i)
      { return Vector<T>::operator[](i); }

    T operator[](size_t i) const
      { return Vector<T>::operator[](i); }

    T& operator()(size_t i)
      { return Vector<T>::operator[](i); }

    T operator()(size_t i) const
      { return Vector<T>::operator[](i); }


    T& operator()(size_t i1, size_t i2, size_t i3, size_t i4)
      {
        if(qn_constrained == 0){
          size_t d1 = _b1.dim();
          size_t d2 = _b2.dim();
          size_t d3 = _b3.dim();
          size_t d4 = _b4.dim();

          return Vector<T>::operator[](i1*d2*d3*d4+i2*d3*d4+i3*d4+i4);
        }else
          return ref_constrained(i1,i2,i3,i4);
      }

    T operator()(size_t i1, size_t i2, size_t i3, size_t i4) const
      {
        if(qn_constrained == 0){
          size_t d1 = _b1.dim();
          size_t d2 = _b2.dim();
          size_t d3 = _b3.dim();
          size_t d4 = _b4.dim();

          return Vector<T>::operator[](i1*d2*d3*d4+i2*d3*d4+i3*d4+i4);
        }else
          return ref_constrained(i1,i2,i3,i4);
      }

    T& operator()(size_t i1, size_t i2, size_t i3, size_t i4, int mask)
      {
        if(mask & MASK_BLOCK1 && mask & MASK_BLOCK2)
          return operator()(i1,i2,i3,i4);
        else if(mask & MASK_BLOCK2 && mask & MASK_BLOCK3)
          return operator()(i3,i1,i2,i4);
        else if(mask & MASK_BLOCK3 && mask & MASK_BLOCK4)
          return operator()(i3,i4,i1,i2);
      }

    T operator()(size_t i1, size_t i2, size_t i3, size_t i4, int mask) const
      {
        if(mask & MASK_BLOCK1 && mask & MASK_BLOCK2)
          return operator()(i1,i2,i3,i4);
        else if(mask & MASK_BLOCK2 && mask & MASK_BLOCK3)
          return operator()(i3,i1,i2,i4);
        else if(mask & MASK_BLOCK3 && mask & MASK_BLOCK4)
          return operator()(i3,i4,i1,i2);
      }


    gslice_iter<T> operator()(size_t i1, size_t i2, Range r3, Range r4); 
    cgslice_iter<T> operator()(size_t i1, size_t i2, Range r3, Range r4) const; 
    gslice_iter<T> operator()(size_t i1, Range r2, size_t i3, Range r4); 
    cgslice_iter<T> operator()(size_t i1, Range r2, size_t i3, Range r4) const; 
    gslice_iter<T> operator()(Range r1, size_t i2, size_t i3, Range r4); 
    cgslice_iter<T> operator()(Range r1, size_t i2, size_t i3, Range r4) const; 

    gslice_iter<T> operator()(size_t i1, Range r2, Range r3, size_t i4); 
    cgslice_iter<T> operator()(size_t i1, Range r2, Range r3, size_t i4) const; 
    gslice_iter<T> operator()(Range r1, size_t i2, Range r3, size_t i4); 
    cgslice_iter<T> operator()(Range r1, size_t i2, Range r3, size_t i4) const; 
    gslice_iter<T> operator()(Range r1, Range r2, size_t i3, size_t i4); 
    cgslice_iter<T> operator()(Range r1, Range r2, size_t i3, size_t i4) const; 

    gslice_iter<T> operator()(Range r1, Range r2, size_t i1, size_t i2, int mask); 
    cgslice_iter<T> operator()(Range r1, Range r2, size_t i1, size_t i2, int mask) const; 

    gslice_iter<T> operator()(size_t i1, size_t i2, QN qn3, QN qn4); 
    cgslice_iter<T> operator()(size_t i1, size_t i2, QN qn3, QN qn4) const; 
    gslice_iter<T> operator()(size_t i1, QN qn2, QN qn3, size_t i4); 
    cgslice_iter<T> operator()(size_t i1, QN qn2, QN qn3, size_t i4) const; 
    gslice_iter<T> operator()(QN qn1, QN qn2, size_t i3, size_t i4); 
    cgslice_iter<T> operator()(QN qn1, QN qn2, size_t i3, size_t i4) const; 

    gslice_iter<T> operator()(QN qn1, QN qn2, size_t i1, size_t i2, int mask); 
    cgslice_iter<T> operator()(QN qn1, QN qn2, size_t i1, size_t i2, int mask) const; 

    slice_iter<T> operator()(size_t i1, size_t i2, size_t i3, Range r4); 
    cslice_iter<T> operator()(size_t i1, size_t i2, size_t i3, Range r4) const; 
    slice_iter<T> operator()(size_t i1, size_t i2, Range r3, size_t i4); 
    cslice_iter<T> operator()(size_t i1, size_t i2, Range r3, size_t i4) const; 
    slice_iter<T> operator()(size_t i1, Range r2, size_t i3, size_t i4); 
    cslice_iter<T> operator()(size_t i1, Range r2, size_t i3, size_t i4) const; 
    slice_iter<T> operator()(Range r1, size_t i2, size_t i3, size_t i4); 
    cslice_iter<T> operator()(Range r1, size_t i2, size_t i3, size_t i4) const; 

    slice_iter<T> operator()(size_t i1, size_t i2, size_t i3, QN qn4); 
    cslice_iter<T> operator()(size_t i1, size_t i2, size_t i3, QN qn4) const; 
    slice_iter<T> operator()(size_t i1, size_t i2, QN qn3, size_t i4); 
    cslice_iter<T> operator()(size_t i1, size_t i2, QN qn3, size_t i4) const; 
    slice_iter<T> operator()(size_t i1, QN qn2, size_t i3, size_t i4); 
    cslice_iter<T> operator()(size_t i1, QN qn2, size_t i3, size_t i4) const; 
    slice_iter<T> operator()(QN qn1, size_t i2, size_t i3, size_t i4); 
    cslice_iter<T> operator()(QN qn1, size_t i2, size_t i3, size_t i4) const; 

    slice_iter<T> operator()(QN qn, size_t i1, size_t i2, size_t i3, int mask); 
    cslice_iter<T> operator()(QN qn, size_t i1, size_t i2, size_t i3, int mask) const; 

    state_slice<T> operator()(QN qn1, QN qn2, QN qn3, QN qn4);
    cstate_slice<T> operator()(QN qn1, QN qn2, QN qn3, QN qn4) const;
    state_slice<T> operator()(const StateSpace &s);
    cstate_slice<T> operator()(const StateSpace &s) const;

    VectorState& resize(const Basis &b1, const Basis &b2, 
                        const Basis &b3, const Basis &b4)
      { 
         Basis aux1(b1), aux2(b2), aux3(b3), aux4(b4);
         init(aux1.pack(),aux2.pack(),aux3.pack(),aux4.pack());
         resize_constrained();
         return *this;
      } 

    VectorState& resize(size_t d1, size_t d2, size_t d3, size_t d4) 
      { 
         init(Basis(d1).pack(),Basis(d2).pack(),Basis(d3).pack(),Basis(d4).pack());
         resize_constrained();
         return *this;
      } 

    VectorState& resize(const PackedBasis &b1, const PackedBasis &b2, 
                        const PackedBasis &b3, const PackedBasis &b4)
      { 
         init(b1,b2,b3,b4);
         resize_constrained();
         return *this;
      } 

    VectorState &constrain(const QN& qn, bool grand_canonical = false)
      {
         _qn = qn;
         qn_constrained = grand_canonical ? 0 : QN::default_mask();
         return *this;
      }

    VectorState &set_qn_mask(const QN& qn, int qn_mask = QN::default_mask())
      {
         _qn = qn;
         qn_constrained = qn_mask;
         return *this;
      }

    VectorState &resize(int qn_mask = QN::default_mask()) { qn_constrained = qn_mask; resize_constrained(); return *this; }
 
    const PackedBasis& b1() const { return _b1; }
    const PackedBasis& b2() const { return _b2; }
    const PackedBasis& b3() const { return _b3; }
    const PackedBasis& b4() const { return _b4; }

    iterator subspace_begin() { return qn_space.begin(); }
    const_iterator subspace_begin() const { return qn_space.begin(); }
    iterator subspace_end() { return qn_space.end(); }
    const_iterator subspace_end() const { return qn_space.end(); }

    bool constrained() const { return (qn_constrained != 0); }
    int qn_mask() const { return qn_constrained; }
    QN qn() const { return _qn; }

    BMatrix<T> density_matrix(int position) const;

    StateSpace get_qn_space(QN qn1, QN qn2, QN qn3, QN qn4) const;
    StateSpace get_qn_space(size_t i1, size_t i2, size_t i3, size_t i4) const;

    VectorState condense(int mask) const;
    VectorState decondense(int mask, const VectorState<T> &orig) const;

//  Streams

    void write(std::ostream &s) const
    {
      _b1.write(s);
      _b2.write(s);
      _b3.write(s);
      _b4.write(s);
      _qn.write(s);
      s.write((const char *)&qn_constrained, sizeof(int));
      Vector<T>::write(s);
    }

    void read(std::istream &s) 
    {
      _b1.read(s);
      _b2.read(s);
      _b3.read(s);
      _b4.read(s);
      _qn.read(s);
      s.read((char *)&qn_constrained, sizeof(int));
      init(_b1,_b2,_b3,_b4);
      if(qn_constrained) resize_constrained();
      Vector<T>::read(s);
    }

};

template<class T>
void
VectorState<T>::resize_constrained()
{
  PackedBasis::const_iterator iter1;          
  PackedBasis::const_iterator iter2;          
  PackedBasis::const_iterator iter3;          
  PackedBasis::const_iterator iter4;          
  size_t dim = 0;

  qn_space.clear();

  for(iter1 = _b1.begin(); iter1 != _b1.end(); iter1++){
    for(iter2 = _b2.begin(); iter2 != _b2.end(); iter2++){
      for(iter3 = _b3.begin(); iter3 != _b3.end(); iter3++){
        for(iter4 = _b4.begin(); iter4 != _b4.end(); iter4++){

          QN qn = (*iter1).qn()+(*iter2).qn()+(*iter3).qn()+(*iter4).qn();

          if(qn_constrained  == 0 || (qn_constrained && qn.equal(_qn,qn_constrained))){
            StateSpace s(*iter1,*iter2,*iter3,*iter4,dim);
            dim += s.dim();

            qn_space.push_back(s);
          }
        }
      }
    }
  }

  Vector<T>::resize(dim);
}

template<class T>
void
VectorState<T>::copy_constrained(const Vector<T>& v)
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  const_iterator iter;
  for(iter = subspace_begin(); iter != subspace_end(); iter++){
    const StateSpace &s = *iter;
    const SubSpace &s1 = s[1];
    const SubSpace &s2 = s[2];
    const SubSpace &s3 = s[3];
    const SubSpace &s4 = s[4];
    for(int i1 = s1.begin(); i1 <= s1.end(); i1++)
      for(int i2 = s2.begin(); i2 <= s2.end(); i2++)
        for(int i3 = s3.begin(); i3 <= s3.end(); i3++)
          for(int i4 = s4.begin(); i4 <= s4.end(); i4++){
            int i = i1*d2*d3*d4+i2*d3*d4+i3*d4+i4;
            operator[](i) = v(i);
          }
  }
}

template<class T>
T&
VectorState<T>::ref_constrained(size_t i1, size_t i2, size_t i3, size_t i4)
{
/*
  QN qn1 = _b1(i1).qn();
  QN qn2 = _b2(i2).qn();
  QN qn3 = _b3(i3).qn();
  QN qn4 = _b4(i4).qn();
*/
  StateSpace s = get_qn_space(i1,i2,i3,i4);
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  size_t index = (i1-s1.begin())*d2*d3*d4+
                 (i2-s2.begin())*d3*d4+
                 (i3-s3.begin())*d4+
                 (i4-s4.begin());
  return Vector<T>::operator[](s.start() + index);
}

template<class T>
T
VectorState<T>::ref_constrained(size_t i1, size_t i2, size_t i3, size_t i4)const
{
/*
  QN qn1 = _b1(i1).qn();
  QN qn2 = _b2(i2).qn();
  QN qn3 = _b3(i3).qn();
  QN qn4 = _b4(i4).qn();
*/
  StateSpace s = get_qn_space(i1,i2,i3,i4);
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  size_t index = (i1-s1.begin())*d2*d3*d4+
                 (i2-s2.begin())*d3*d4+
                 (i3-s3.begin())*d4+
                 (i4-s4.begin());
  return Vector<T>::operator[](s.start() + index);
}

/////////////////////////////////////////////////////////////////////////
template<class T>
gslice_iter<T> 
VectorState<T>::operator()(size_t i1, size_t i2, Range r3, Range r4)
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s3(i1*d2*d3*d4+i2*d3*d4+r3.start()*d4,r3.size(),r3.stride()*d4);
  slice s4(r4.start(),r4.size(),r4.stride());
  return gslice_iter<T>(this, s3, s4);
}

template<class T>
gslice_iter<T> 
VectorState<T>::operator()(size_t i1, Range r2, Range r3, size_t i4)
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s2(i1*d2*d3*d4+r2.start()*d3*d4,r2.size(),r2.stride()*d3*d4);
  slice s3(r3.start()*d4+i4,r3.size(),r3.stride()*d4);
  return gslice_iter<T>(this, s2, s3);
}

template<class T>
gslice_iter<T> 
VectorState<T>::operator()(Range r1, Range r2, size_t i3, size_t i4)
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s1(r1.start()*d2*d3*d4,r1.size(),r1.stride()*d2*d3*d4);
  slice s2(r2.start()*d3*d4+i3*d4+i4,r2.size(),r2.stride()*d3*d4);
  return gslice_iter<T>(this, s1, s2);
}

template<class T>
gslice_iter<T> 
VectorState<T>::operator()(size_t i1, Range r2, size_t i3, Range r4)
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s2(i1*d2*d3*d4+r2.start()*d3*d4+i3*d4,r2.size(),r2.stride()*d3*d4);
  slice s4(r4.start(),r4.size(),r4.stride());
  return gslice_iter<T>(this, s2, s4);
}

template<class T>
gslice_iter<T> 
VectorState<T>::operator()(Range r1, size_t i2, size_t i3, Range r4)
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s1(r1.start()*d2*d3*d4+i2*d3*d4+i3*d4,r1.size(),r1.stride()*d2*d3*d4);
  slice s4(r4.start(),r4.size(),r4.stride());
  return gslice_iter<T>(this, s1, s4);
}

template<class T>
gslice_iter<T> 
VectorState<T>::operator()(Range r1, size_t i2, Range r3, size_t i4)
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s1(r1.start()*d2*d3*d4+i2*d3*d4,r1.size(),r1.stride()*d2*d3*d4);
  slice s3(r3.start()*d4+i4,r3.size(),r3.stride()*d4);
  return gslice_iter<T>(this, s1, s3);
}

///////////////////////////////////////////////////////////////////////////

template<class T>
cgslice_iter<T> 
VectorState<T>::operator()(size_t i1, size_t i2, Range r3, Range r4) const
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s3(i1*d2*d3*d4+i2*d3*d4+r3.start()*d4,r3.size(),r3.stride()*d4);
  slice s4(r4.start(),r4.size(),r4.stride());
  return cgslice_iter<T>(this, s3, s4);
}

template<class T>
cgslice_iter<T> 
VectorState<T>::operator()(size_t i1, Range r2, Range r3, size_t i4) const
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s2(i1*d2*d3*d4+r2.start()*d3*d4,r2.size(),r2.stride()*d3*d4);
  slice s3(r3.start()*d4+i4,r3.size(),r3.stride()*d4);
  return cgslice_iter<T>(this, s2, s3);
}

template<class T>
cgslice_iter<T> 
VectorState<T>::operator()(size_t i1, Range r2, size_t i3, Range r4) const
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s2(i1*d2*d3*d4+r2.start()*d3*d4+i3*d4,r2.size(),r2.stride()*d3*d4);
  slice s4(r4.start(),r4.size(),r4.stride());
  return cgslice_iter<T>(this, s2, s4);
}

template<class T>
cgslice_iter<T> 
VectorState<T>::operator()(Range r1, Range r2, size_t i3, size_t i4) const
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s1(r1.start()*d2*d3*d4,r1.size(),r1.stride()*d2*d3*d4);
  slice s2(r2.start()*d3*d4+i3*d4+i4,r2.size(),r2.stride()*d3*d4);
  return cgslice_iter<T>(this, s1, s2);
}

template<class T>
cgslice_iter<T> 
VectorState<T>::operator()(Range r1, size_t i2, size_t i3, Range r4) const
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s1(r1.start()*d2*d3*d4+i2*d3*d4+i3*d4,r1.size(),r1.stride()*d2*d3*d4);
  slice s4(r4.start(),r4.size(),r4.stride());
  return cgslice_iter<T>(this, s1, s4);
}

template<class T>
cgslice_iter<T> 
VectorState<T>::operator()(Range r1, size_t i2, Range r3, size_t i4) const
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s1(r1.start()*d2*d3*d4+i2*d3*d4,r1.size(),r1.stride()*d2*d3*d4);
  slice s3(r3.start()*d4+i4,r3.size(),r3.stride()*d4);
  return cgslice_iter<T>(this, s1, s3);
}

////////////////////////////////////////////////////////
template<class T>
slice_iter<T> 
VectorState<T>::operator()(size_t i1, size_t i2, size_t i3, Range r4)
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s4(i1*d2*d3*d4+i2*d3*d4+i3*d4+r4.start(),r4.size(),r4.stride());
  return slice_iter<T>(this, s4);
}

template<class T>
slice_iter<T> 
VectorState<T>::operator()(size_t i1, size_t i2, Range r3, size_t i4)
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s3(i1*d2*d3*d4+i2*d3*d4+r3.start()*d4+i4,r3.size(),r3.stride()*d4);
  return slice_iter<T>(this, s3);
}

template<class T>
slice_iter<T> 
VectorState<T>::operator()(size_t i1, Range r2, size_t i3, size_t i4)
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s2(i1*d2*d3*d4+r2.start()*d3*d4+i3*d4+i4,r2.size(),r2.stride()*d3*d4);
  return slice_iter<T>(this, s2);
}

template<class T>
slice_iter<T> 
VectorState<T>::operator()(Range r1, size_t i2, size_t i3, size_t i4)
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s1(r1.start()*d2*d3*d4+i2*d3*d4+i3*d4+i4,r1.size(),r1.stride()*d2*d3*d4);
  return slice_iter<T>(this, s1);
}

template<class T>
cslice_iter<T> 
VectorState<T>::operator()(size_t i1, size_t i2, size_t i3, Range r4) const
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s4(i1*d2*d3*d4+i2*d3*d4+i3*d4+r4.start(),r4.size(),r4.stride());
  return cslice_iter<T>(this, s4);
}

template<class T>
cslice_iter<T> 
VectorState<T>::operator()(size_t i1, size_t i2, Range r3, size_t i4) const
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s3(i1*d2*d3*d4+i2*d3*d4+r3.start()*d4+i4,r3.size(),r3.stride()*d4);
  return cslice_iter<T>(this, s3);
}

template<class T>
cslice_iter<T> 
VectorState<T>::operator()(size_t i1, Range r2, size_t i3, size_t i4) const
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s2(i1*d2*d3*d4+r2.start()*d3*d4+i3*d4+i4,r2.size(),r2.stride()*d3*d4);
  return cslice_iter<T>(this, s2);
}

template<class T>
cslice_iter<T> 
VectorState<T>::operator()(Range r1, size_t i2, size_t i3, size_t i4) const
{
  size_t d1 = _b1.dim();
  size_t d2 = _b2.dim();
  size_t d3 = _b3.dim();
  size_t d4 = _b4.dim();

  if(qn_constrained != 0) cerr << "*** WARNING: Constrained StateVector. This operation may not work properly\n";

  slice s1(r1.start()*d2*d3*d4+i2*d3*d4+i3*d4+i4,r1.size(),r1.stride()*d2*d3*d4);
  return cslice_iter<T>(this, s1);
}

////////////////////////////////////////////////////////

template<class T>
gslice_iter<T> 
VectorState<T>::operator()(Range r1, Range r2, size_t i1, size_t i2, int mask) 
{
  if((mask & MASK_BLOCK1) && (mask & MASK_BLOCK2))
    return operator()(r1, r2, i1, i2); 
  else if((mask & MASK_BLOCK2) && (mask & MASK_BLOCK3))
    return operator()(i1, r1, r2, i2); 
  else if((mask & MASK_BLOCK3) && (mask & MASK_BLOCK4))
    return operator()(i1, i2, r1, r2); 
  else if((mask & MASK_BLOCK2) && (mask & MASK_BLOCK4))
    return operator()(i1, r1, i2, r2); 
  else if((mask & MASK_BLOCK1) && (mask & MASK_BLOCK4))
    return operator()(r1, i1, i2, r2); 
  else if((mask & MASK_BLOCK1) && (mask & MASK_BLOCK3))
    return operator()(r1, i1, r2, i2); 

  return gslice_iter<T>();
}

template<class T>
cgslice_iter<T> 
VectorState<T>::operator()(Range r1, Range r2, size_t i1, size_t i2, int mask) const
{
  if((mask & MASK_BLOCK1) && (mask & MASK_BLOCK2))
    return operator()(r1, r2, i1, i2); 
  else if((mask & MASK_BLOCK2) && (mask & MASK_BLOCK3))
    return operator()(i1, r1, r2, i2); 
  else if((mask & MASK_BLOCK3) && (mask & MASK_BLOCK4))
    return operator()(i1, i2, r1, r2); 
  else if((mask & MASK_BLOCK2) && (mask & MASK_BLOCK4))
    return operator()(i1, r1, i2, r2); 
  else if((mask & MASK_BLOCK1) && (mask & MASK_BLOCK4))
    return operator()(r1, i1, i2, r2); 
  else if((mask & MASK_BLOCK1) && (mask & MASK_BLOCK3))
    return operator()(r1, i1, r2, i2); 

  return cgslice_iter<T>();
}

/////////////////////////////////////////////////////////////////////////
template<class T>
StateSpace 
VectorState<T>::get_qn_space(QN qn1, QN qn2, QN qn3, QN qn4) const
{
/*
  if(qn_constrained == 0){
    typename VectorState<T>::const_iterator iter;
    for(iter = subspace_begin(); iter != subspace_end(); iter++){
      const SubSpace &s1 = (*iter)[1];
      const SubSpace &s2 = (*iter)[2];
      const SubSpace &s3 = (*iter)[3];
      const SubSpace &s4 = (*iter)[4];
 
      if(s1.qn() == qn1 && s2.qn() == qn2 && s3.qn() == qn3 && s4.qn() == qn4)
        return (*iter);
    }
//    cout << "*** ERROR 1: Indices out of constraint\n";
    StateSpace error;
    error._start = -1;
    return error;
  }
*/

  int origin = 0;
  int end = qn_space.size() - 1;
  int index = -1;

  while(origin <= end){
    int index_old = index;
    index = (origin + end) / 2;
 
    const StateSpace &s = qn_space[index];

    QN qn = s[1].qn() + s[2].qn() + s[3].qn() + s[4].qn();

    if(s[1].qn() == qn1 && s[2].qn() == qn2 && 
                           s[3].qn() == qn3 && s[4].qn() == qn4) return s; 
    if(s[1].qn() > qn1 ||
       (s[1].qn() == qn1 && s[2].qn() > qn2) ||
       (s[1].qn() == qn1 && s[2].qn() == qn2 && s[3].qn() > qn3) ||
       (s[1].qn() == qn1 && s[2].qn() == qn2 && s[3].qn() == qn3 && s[4].qn() > qn4)) 
      end = index;
    else
      origin = index;

    if(index == index_old){
      if(index == end)
        end = end - 1;
      else
        origin = origin + 1;
    } 
 
  } 
//  cout << "*** ERROR 2: Indices out of constraint\n";

  StateSpace error;
  error._start = -1;
  return error;
}

template<class T>
StateSpace 
VectorState<T>::get_qn_space(size_t i1, size_t i2, size_t i3, size_t i4) const
{
  QN qn1, qn2, qn3, qn4;
  qn1 = _b1(i1).qn();
  qn2 = _b2(i2).qn();
  qn3 = _b3(i3).qn();
  qn4 = _b4(i4).qn();
  return(get_qn_space(qn1,qn2,qn3,qn4));

//-------------------------
  const_iterator iter;
  for(iter = subspace_begin(); iter != subspace_end(); iter++){
    const SubSpace &s1 = (*iter)[1];
    const SubSpace &s2 = (*iter)[2];
    const SubSpace &s3 = (*iter)[3];
    const SubSpace &s4 = (*iter)[4];

    if(s1.begin() <= i1 && i1 <= s1.end() &&
       s2.begin() <= i2 && i2 <= s2.end() &&
       s3.begin() <= i3 && i3 <= s3.end() &&
       s4.begin() <= i4 && i4 <= s4.end())
                  return (*iter);
  }
//  cout << "*** ERROR 3: Indices out of constraint\n";
  StateSpace error;
  error._start = -1;
  return error;
}

/////////////////////////////////////////////////////////////////////////

template<class T>
gslice_iter<T> 
VectorState<T>::operator()(size_t i1, size_t i2, QN qn3, QN qn4)
{

  if(!qn_constrained){
    SubSpace s3 = _b3(qn3);
    SubSpace s4 = _b4(qn4);
    return operator()(i1,i2,s3,s4);
  }
    
  QN qn1 = _b1(i1).qn();
  QN qn2 = _b2(i2).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4); 
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t start = (i1 - s1.begin())*d2*s3.dim()*s4.dim() +
                 (i2 - s2.begin())*s3.dim()*s4.dim();
  slice _s3(s.start()+start,s3.dim(),s4.dim());
  slice _s4(0,s4.dim(),1);
  return gslice_iter<T>(this, _s3, _s4);
}

template<class T>
gslice_iter<T> 
VectorState<T>::operator()(size_t i1, QN qn2, QN qn3, size_t i4)
{

  if(!qn_constrained){
    SubSpace s2 = _b2(qn2);
    SubSpace s3 = _b3(qn3);
    return operator()(i1,s2,s3,i4);
  }
    
  QN qn1 = _b1(i1).qn();
  QN qn4 = _b4(i4).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4); 
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d4 = s4.dim();
  size_t start = (i1 - s1.begin())*s2.dim()*s3.dim()*d4 + (i4 - s4.begin());
  slice _s2(s.start()+start,s2.dim(),s3.dim()*d4);
  slice _s3(0,s3.dim(),d4);
  return gslice_iter<T>(this, _s2, _s3);
}

template<class T>
gslice_iter<T> 
VectorState<T>::operator()(QN qn1, QN qn2, size_t i3, size_t i4)
{

  if(!qn_constrained){
    SubSpace s1 = _b1(qn1);
    SubSpace s2 = _b2(qn2);
    return operator()(s1,s2,i3,i4);
  }
    
  QN qn3 = _b3(i3).qn();
  QN qn4 = _b4(i4).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4); 
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  size_t start = (i3 - s3.begin())*d4 + (i4 - s4.begin());
  slice _s1(s.start()+start,s1.dim(),s2.dim()*d3*d4);
  slice _s2(0,s2.dim(),d3*d4);
  return gslice_iter<T>(this, _s1, _s2);
}

template<class T>
cgslice_iter<T> 
VectorState<T>::operator()(size_t i1, size_t i2, QN qn3, QN qn4) const
{

  if(!qn_constrained){
    SubSpace s3 = _b3(qn3);
    SubSpace s4 = _b4(qn4);
    return operator()(i1,i2,s3,s4);
  }
    
  QN qn1 = _b1(i1).qn();
  QN qn2 = _b2(i2).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4); 
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t start = (i1 - s1.begin())*d2*s3.dim()*s4.dim() +
                 (i2 - s2.begin())*s3.dim()*s4.dim();
  slice _s3(s.start()+start,s3.dim(),s4.dim());
  slice _s4(0,s4.dim(),1);
  return cgslice_iter<T>(this, _s3, _s4);
}

template<class T>
cgslice_iter<T> 
VectorState<T>::operator()(size_t i1, QN qn2, QN qn3, size_t i4) const
{

  if(!qn_constrained){
    SubSpace s2 = _b2(qn2);
    SubSpace s3 = _b3(qn3);
    return operator()(i1,s2,s3,i4);
  }
    
  QN qn1 = _b1(i1).qn();
  QN qn4 = _b4(i4).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4); 
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d4 = s4.dim();
  size_t start = (i1 - s1.begin())*s2.dim()*s3.dim()*d4 + (i4 - s4.begin());
  slice _s2(s.start()+start,s2.dim(),s3.dim()*d4);
  slice _s3(0,s3.dim(),d4);
  return cgslice_iter<T>(this, _s2, _s3);
}

template<class T>
cgslice_iter<T> 
VectorState<T>::operator()(QN qn1, QN qn2, size_t i3, size_t i4) const
{

  if(!qn_constrained){
    SubSpace s1 = _b1(qn1);
    SubSpace s2 = _b2(qn2);
    return operator()(s1,s2,i3,i4);
  }
    
  QN qn3 = _b3(i3).qn();
  QN qn4 = _b4(i4).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4); 
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  size_t start = (i3 - s3.begin())*d4 + (i4 - s4.begin());
  slice _s1(s.start()+start,s1.dim(),s2.dim()*d3*d4);
  slice _s2(0,s2.dim(),d3*d4);
  return cgslice_iter<T>(this, _s1, _s2);
}

////////////////////////////////////////////////////////

template<class T>
gslice_iter<T> 
VectorState<T>::operator()(QN qn1, QN qn2, size_t i1, size_t i2, int mask) 
{
  if((mask & MASK_BLOCK1) && (mask & MASK_BLOCK2))
    return operator()(qn1, qn2, i1, i2); 
  else if((mask & MASK_BLOCK2) && (mask & MASK_BLOCK3))
    return operator()(i1, qn1, qn2, i2); 
  else if((mask & MASK_BLOCK3) && (mask & MASK_BLOCK4))
    return operator()(i1, i2, qn1, qn2); 

  return gslice_iter<T>();
}

template<class T>
cgslice_iter<T> 
VectorState<T>::operator()(QN qn1, QN qn2, size_t i1, size_t i2, int mask) const
{
  if((mask & MASK_BLOCK1) && (mask & MASK_BLOCK2))
    return operator()(qn1, qn2, i1, i2); 
  else if((mask & MASK_BLOCK2) && (mask & MASK_BLOCK3))
    return operator()(i1, qn1, qn2, i2); 
  else if((mask & MASK_BLOCK3) && (mask & MASK_BLOCK4))
    return operator()(i1, i2, qn1, qn2); 

  return cgslice_iter<T>();
}

////////////////////////////////////////////////////////

template<class T>
slice_iter<T> 
VectorState<T>::operator()(QN qn1, size_t i2, size_t i3, size_t i4)
{

  if(!qn_constrained){
    SubSpace s1 = _b1(qn1);
    return operator()(s1,i2,i3,i4);
  }
    
  QN qn2 = _b2(i2).qn();
  QN qn3 = _b3(i3).qn();
  QN qn4 = _b4(i4).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4);
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  size_t start = (i2 - s2.begin())*d3*d4 +
                 (i3 - s3.begin())*d4 + (i4 - s4.begin());
  slice _s1(s.start()+start,d1,d2*d3*d4);
  return slice_iter<T>(this, _s1);
}

template<class T>
slice_iter<T> 
VectorState<T>::operator()(size_t i1, QN qn2, size_t i3, size_t i4)
{

  if(!qn_constrained){
    SubSpace s2 = _b2(qn2);
    return operator()(i1,s2,i3,i4);
  }
    
  QN qn1 = _b1(i1).qn();
  QN qn3 = _b3(i3).qn();
  QN qn4 = _b4(i4).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4);
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  size_t start = (i1 - s1.begin())*d2*d3*d4 +
                 (i3 - s3.begin())*d4 + (i4 - s4.begin());
  slice _s2(s.start()+start,d2,d3*d4);
  return slice_iter<T>(this, _s2);
}

template<class T>
slice_iter<T> 
VectorState<T>::operator()(size_t i1, size_t i2, QN qn3, size_t i4)
{

  if(!qn_constrained){
    SubSpace s3 = _b3(qn3);
    return operator()(i1,i2,s3,i4);
  }
    
  QN qn1 = _b1(i1).qn();
  QN qn2 = _b2(i2).qn();
  QN qn4 = _b4(i4).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4);
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  size_t start = (i1 - s1.begin())*d2*d3*d4 +
                 (i2 - s2.begin())*d3*d4 + 
                 (i4 - s4.begin());
  slice _s3(s.start()+start,d3,d4);
  return slice_iter<T>(this, _s3);
}

template<class T>
slice_iter<T> 
VectorState<T>::operator()(size_t i1, size_t i2, size_t i3, QN qn4)
{
  if(!qn_constrained){
    SubSpace s4 = _b4(qn4);
    return operator()(i1,i2,i3,s4);
  }
    
  QN qn1 = _b1(i1).qn();
  QN qn2 = _b2(i2).qn();
  QN qn3 = _b3(i3).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4);
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  size_t start = (i1 - s1.begin())*d2*d3*d4 +
                 (i2 - s2.begin())*d3*d4 + 
                 (i3 - s3.begin())*d4;
  slice _s4(s.start()+start,d4,1);
  return slice_iter<T>(this, _s4);
}

template<class T>
cslice_iter<T> 
VectorState<T>::operator()(QN qn1, size_t i2, size_t i3, size_t i4) const
{

  if(!qn_constrained){
    SubSpace s1 = _b1(qn1);
    return operator()(s1,i2,i3,i4);
  }
    
  QN qn2 = _b2(i2).qn();
  QN qn3 = _b3(i3).qn();
  QN qn4 = _b4(i4).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4);
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  size_t start = (i2 - s2.begin())*d3*d4 +
                 (i3 - s3.begin())*d4 + (i4 - s4.begin());
  slice _s1(s.start()+start,d1,d2*d3*d4);
  return cslice_iter<T>(this, _s1);
}

template<class T>
cslice_iter<T> 
VectorState<T>::operator()(size_t i1, QN qn2, size_t i3, size_t i4) const
{

  if(!qn_constrained){
    SubSpace s2 = _b2(qn2);
    return operator()(i1,s2,i3,i4);
  }
    
  QN qn1 = _b1(i1).qn();
  QN qn3 = _b3(i3).qn();
  QN qn4 = _b4(i4).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4);
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  size_t start = (i1 - s1.begin())*d2*d3*d4 +
                 (i3 - s3.begin())*d4 + (i4 - s4.begin());
  slice _s2(s.start()+start,d2,d3*d4);
  return cslice_iter<T>(this, _s2);
}

template<class T>
cslice_iter<T> 
VectorState<T>::operator()(size_t i1, size_t i2, QN qn3, size_t i4) const
{

  if(!qn_constrained){
    SubSpace s3 = _b3(qn3);
    return operator()(i1,i2,s3,i4);
  }
    
  QN qn1 = _b1(i1).qn();
  QN qn2 = _b2(i2).qn();
  QN qn4 = _b4(i4).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4);
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  size_t start = (i1 - s1.begin())*d2*d3*d4 +
                 (i2 - s2.begin())*d3*d4 + 
                 (i4 - s4.begin());
  slice _s3(s.start()+start,d3,d4);
  return cslice_iter<T>(this, _s3);
}

template<class T>
cslice_iter<T> 
VectorState<T>::operator()(size_t i1, size_t i2, size_t i3, QN qn4) const
{
  if(!qn_constrained){
    SubSpace s4 = _b4(qn4);
    return operator()(i1,i2,i3,s4);
  }
    
  QN qn1 = _b1(i1).qn();
  QN qn2 = _b2(i2).qn();
  QN qn3 = _b3(i3).qn();
  StateSpace s = get_qn_space(qn1,qn2,qn3,qn4);
  if(s.start() == -1) cout << "ERROR: VectorState\n";
  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  size_t start = (i1 - s1.begin())*d2*d3*d4 +
                 (i2 - s2.begin())*d3*d4 + 
                 (i3 - s3.begin())*d4;
  slice _s4(s.start()+start,d4,1);
  return cslice_iter<T>(this, _s4);
}

////////////////////////////////////////////////////////

template<class T>
slice_iter<T> 
VectorState<T>::operator()(QN qn, size_t i1, size_t i2, size_t i3, int mask) 
{
  if(mask & MASK_BLOCK1)
    return operator()(qn, i1, i2, i3); 
  else if(mask & MASK_BLOCK2)
    return operator()(i1, qn, i2, i3); 
  else if(mask & MASK_BLOCK3)
    return operator()(i1, i2, qn, i3); 
  else if(mask & MASK_BLOCK4)
    return operator()(i1, i2, i3, qn); 

  return slice_iter<T>();
}

template<class T>
cslice_iter<T> 
VectorState<T>::operator()(QN qn, size_t i1, size_t i2, size_t i3, int mask) const 
{
  if(mask & MASK_BLOCK1)
    return operator()(qn, i1, i2, i3); 
  else if(mask & MASK_BLOCK2)
    return operator()(i1, qn, i2, i3); 
  else if(mask & MASK_BLOCK3)
    return operator()(i1, i2, qn, i3); 
  else if(mask & MASK_BLOCK4)
    return operator()(i1, i2, i3, qn); 

  return cslice_iter<T>();
}

////////////////////////////////////////////////////////

template<class T>
state_slice<T> 
VectorState<T>::operator()(const StateSpace &s)
{
  if(s.start() == -1) return state_slice<T>();

  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  slice _s1(s.start(),d1,d2*d3*d4);
  slice _s2(0,d2,d3*d4);
  slice _s3(0,d3,d4);
  slice _s4(0,d4,1);
  return state_slice<T>(this, _s1, _s2, _s3, _s4);
}

template<class T>
cstate_slice<T> 
VectorState<T>::operator()(const StateSpace &s) const
{
  if(s.start() == -1) return cstate_slice<T>();

  const SubSpace &s1 = s[1];
  const SubSpace &s2 = s[2];
  const SubSpace &s3 = s[3];
  const SubSpace &s4 = s[4];
  size_t d1 = s1.dim();
  size_t d2 = s2.dim();
  size_t d3 = s3.dim();
  size_t d4 = s4.dim();
  slice _s1(s.start(),d1,d2*d3*d4);
  slice _s2(0,d2,d3*d4);
  slice _s3(0,d3,d4);
  slice _s4(0,d4,1);
  return cstate_slice<T>(this, _s1, _s2, _s3, _s4);
}

template<class T>
state_slice<T> 
VectorState<T>::operator()(QN qn1, QN qn2, QN qn3, QN qn4)
{
  const StateSpace &s = get_qn_space(qn1,qn2,qn3,qn4);
  return operator()(s);
}

template<class T>
cstate_slice<T> 
VectorState<T>::operator()(QN qn1, QN qn2, QN qn3, QN qn4) const
{
  const StateSpace &s = get_qn_space(qn1,qn2,qn3,qn4);
  return operator()(s);
}

////////////////////////////////////////////////////////

template<class T>
BMatrix<T>
VectorState<T>::density_matrix(int position) const
{
  Basis basis(_b1,_b2);
  if(position == RIGHT) basis = Basis(_b3,_b4);
  basis.reorder();
  BMatrix<T> rho(basis);

  PackedBasis::iterator biter;
  for(biter = basis.subspace_begin(); biter != basis.subspace_end(); biter++){
    SubSpace s = *biter;
    SubMatrix<T> block(s.qn(),s,s);
    rho.push_back(block);
  }

  VectorState aux;
  aux = condense1(position == LEFT ? (MASK_BLOCK1|MASK_BLOCK2) : (MASK_BLOCK3|MASK_BLOCK4));

  QN qn0(-999);
  SubMatrix<T> *_sub;
  typename VectorState<T>::iterator iter;

  for(iter = aux.subspace_begin(); iter != aux.subspace_end(); iter++){
    StateSpace s = *iter;
    if(s[1].qn() != qn0 || !_sub || qn0[0] == -999){
      qn0 = s[1].qn();
      _sub = rho.block(qn0);
    }
    if(!_sub) continue;

    SubMatrix<T> &sub = *_sub;
    state_slice<T> s_slice = aux(s);

    slice r3(0,s[3].size(),1);
    slice r4(0,s[4].size(),1);
    Matrix<T> msubi(r3.size(),r4.size());
    Matrix<T> msubj(r3.size(),r4.size());
    int ntot = r3.size()*r4.size();
    for(int i = 0; i < s[1].dim(); i++){
      gslice_iter<T> subi = s_slice(i,0,r3,r4);
      msubi = subi;

      sub(i,i) += dot_product(ntot,msubi.array(),1,msubi.array(),1);

      for(int j = i+1; j < s[1].dim(); j++){ // rho is Hermitian
        gslice_iter<T> subj = s_slice(j,0,r3,r4);
        msubj = subj;

        T sum = dot_product(ntot,msubi.array(),1,msubj.array(),1);

        sub(i,j) += sum;
        sub(j,i) += conj(sum);
      }
    }
  }

  return rho;
}

////////////////////////////////////////////////////////////////////
// condense1:
// this function builds a new state vector merging two blocks.
// The new resulting block is on the left 
// We use this vector for the construction of the density matrix.
////////////////////////////////////////////////////////////////////
template<class T>
VectorState<T>
VectorState<T>::condense1 (int mask) const
{
// we build an auxiliary array of subspaces with the two blocks added.
// this will be used to reorder the indices of the subspaces, considering
// the two blocks as a single block

  Vector<StateSpace> aux_qn_space(qn_space.size());
  typename std::vector<StateSpace>::const_iterator siter;
  typename Vector<StateSpace>::iterator titer;
  for(siter = qn_space.begin(), titer = aux_qn_space.begin(); siter != qn_space.end(); siter++, titer++){
    const StateSpace &s = *siter;
    StateSpace ss;
    if(mask == (MASK_BLOCK1|MASK_BLOCK2))
      ss = StateSpace(SubSpace(s[1].qn()+s[2].qn(),0,0),SubSpace(),s[3],s[4]);
    else if(mask == (MASK_BLOCK2|MASK_BLOCK3))
      ss = StateSpace(SubSpace(s[2].qn()+s[3].qn(),0,0),SubSpace(),s[1],s[4]);
    else if(mask == (MASK_BLOCK3|MASK_BLOCK4))
      ss = StateSpace(SubSpace(s[3].qn()+s[4].qn(),0,0),SubSpace(),s[1],s[2]);
    (*titer) = ss;
  }

// we reorder the subspaces considering the two blocks as a single block
  
  Vector<size_t> idx(aux_qn_space.size());
  indexx<StateSpace, Vector<StateSpace> >(aux_qn_space.size(), aux_qn_space, idx);

// we build the condensed state vector with the two blocks added,
// rearranging the subspaces

  Basis aux_basis(1);
  aux_basis.reorder();

  VectorState<T> aux(*this);
  if(mask == (MASK_BLOCK1|MASK_BLOCK2)){
     Basis basis = Basis(_b1,_b2);
     basis.reorder();
     aux.resize(basis.subspaces(),aux_basis.subspaces(),_b3,_b4); 
  }else if(mask == (MASK_BLOCK2|MASK_BLOCK3)){
     Basis basis = Basis(_b2,_b3);
     basis.reorder();
     aux.resize(basis.subspaces(),aux_basis.subspaces(),_b1,_b4); 
  }else if(mask == (MASK_BLOCK3|MASK_BLOCK4)){
     Basis basis = Basis(_b3,_b4);
     basis.reorder();
     aux.resize(basis.subspaces(),aux_basis.subspaces(),_b1,_b2); 
  }

  typename Vector<T>::iterator aux_iter = aux.begin();
  for(int i = 0; i < idx.size(); i++){
    StateSpace s = qn_space[idx[i]];

    cstate_slice<T> s_slice = operator()(s);
    if(mask == (MASK_BLOCK1|MASK_BLOCK2)){
      for(int i1 = 0; i1 < s[1].dim(); i1++){ 
        for(int i2 = 0; i2 < s[2].dim(); i2++){ 
          for(int i3 = 0; i3 < s[3].dim(); i3++){ 
            for(int i4 = 0; i4 < s[4].dim(); i4++){ 
              *aux_iter = s_slice(i1,i2,i3,i4); 
              aux_iter++;
            }
          }
        }
      }
    }else if(mask == (MASK_BLOCK2|MASK_BLOCK3)){
      for(int i2 = 0; i2 < s[2].dim(); i2++){ 
        for(int i3 = 0; i3 < s[3].dim(); i3++){ 
          for(int i1 = 0; i1 < s[1].dim(); i1++){ 
            for(int i4 = 0; i4 < s[4].dim(); i4++){ 
              *aux_iter = s_slice(i1,i2,i3,i4); 
              aux_iter++;
            }
          }
        }
      }
    }else if(mask == (MASK_BLOCK3|MASK_BLOCK4)){
      for(int i3 = 0; i3 < s[3].dim(); i3++){ 
        for(int i4 = 0; i4 < s[4].dim(); i4++){ 
          for(int i1 = 0; i1 < s[1].dim(); i1++){ 
            for(int i2 = 0; i2 < s[2].dim(); i2++){ 
              *aux_iter = s_slice(i1,i2,i3,i4); 
              aux_iter++;
            }
          }
        }
      }
    }
  }

  return aux;
}

template<class T>
VectorState<T>
VectorState<T>::decondense1 (int mask, const VectorState<T> &orig) const
{
  VectorState<T> aux(orig);

// we build an auxiliary array of subspaces with the two blocks added.
// this will be used to reorder the indices of the subspaces, considering
// the two blocks as a single block

  Vector<StateSpace> aux_qn_space(aux.qn_space.size());
  typename vector<StateSpace>::const_iterator siter;
  typename Vector<StateSpace>::iterator titer;
  for(siter = aux.qn_space.begin(), titer = aux_qn_space.begin(); siter != aux.qn_space.end(); siter++, titer++){
    const StateSpace &s = *siter;
    StateSpace ss;
    if(mask == (MASK_BLOCK1|MASK_BLOCK2))
      ss = StateSpace(SubSpace(s[1].qn()+s[2].qn(),0,0),SubSpace(),s[3],s[4]);
    else if(mask == (MASK_BLOCK2|MASK_BLOCK3))
      ss = StateSpace(SubSpace(s[2].qn()+s[3].qn(),0,0),SubSpace(),s[1],s[4]);
    else if(mask == (MASK_BLOCK3|MASK_BLOCK4))
      ss = StateSpace(SubSpace(s[3].qn()+s[4].qn(),0,0),SubSpace(),s[1],s[2]);
    (*titer) = ss;
  }

// we reorder the subspaces considering the two blocks as a single block
  
  Vector<size_t> idx(aux_qn_space.size());
  indexx<StateSpace, Vector<StateSpace> >(aux_qn_space.size(), aux_qn_space, idx);

// we build the condensed state vector with the two blocks added,
// rearranging the subspaces

  typename Vector<T>::const_iterator aux_iter = Vector<T>::begin();
  for(int i = 0; i < idx.size(); i++){
    StateSpace s = aux.qn_space[idx[i]];

    state_slice<T> s_slice = aux(s);
    if(mask == (MASK_BLOCK1|MASK_BLOCK2)){
      for(int i1 = 0; i1 < s[1].dim(); i1++){ 
        for(int i2 = 0; i2 < s[2].dim(); i2++){ 
          for(int i3 = 0; i3 < s[3].dim(); i3++){ 
            for(int i4 = 0; i4 < s[4].dim(); i4++){ 
              s_slice(i1,i2,i3,i4) = *aux_iter; 
              aux_iter++;
            }
          }
        }
      }
    }else if(mask == (MASK_BLOCK2|MASK_BLOCK3)){
      for(int i2 = 0; i2 < s[2].dim(); i2++){ 
        for(int i3 = 0; i3 < s[3].dim(); i3++){ 
          for(int i1 = 0; i1 < s[1].dim(); i1++){ 
            for(int i4 = 0; i4 < s[4].dim(); i4++){ 
              s_slice(i1,i2,i3,i4) = *aux_iter; 
              aux_iter++;
            }
          }
        }
      }
    }else if(mask == (MASK_BLOCK3|MASK_BLOCK4)){
      for(int i3 = 0; i3 < s[3].dim(); i3++){ 
        for(int i4 = 0; i4 < s[4].dim(); i4++){ 
          for(int i1 = 0; i1 < s[1].dim(); i1++){ 
            for(int i2 = 0; i2 < s[2].dim(); i2++){ 
              s_slice(i1,i2,i3,i4) = *aux_iter; 
              aux_iter++;
            }
          }
        }
      }
    }
  }

  return aux;
}

////////////////////////////////////////////////////////////////////
// condense:
// this function builds a new state vector merging two blocks.
////////////////////////////////////////////////////////////////////
template<class T>
VectorState<T>
VectorState<T>::condense (int mask) const
{
  typename VectorState<T>::const_iterator siter;

  VectorState<T> aux(*this);
  aux.index.resize(aux.size());

  if(mask == (MASK_BLOCK1|MASK_BLOCK2)){
    Basis basis(_b1,_b2);
    basis.reorder();

    Basis aux_basis(1);
    aux_basis.reorder();

    aux.resize(basis.subspaces(),aux_basis.subspaces(),_b3,_b4);

    for(siter = aux.subspace_begin(); siter != aux.subspace_end(); siter++){
      StateSpace ss = *siter;

      state_slice<T> aux_slice = aux(ss);
      for(int i1 = 0; i1 < ss[1].dim(); i1++){
        State sj = basis(i1+ss[1].begin());
        int j1 = sj.i1;
        int j2 = sj.i2;
        StateSpace cs = get_qn_space(sj.qn1(),sj.qn2(),ss[3].qn(),ss[4].qn());
        cstate_slice<T> c_slice = operator()(cs);
 
        for(int i3 = 0; i3 < ss[3].dim(); i3++){
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            size_t index_dest = aux_slice.index(i1,0,i3,i4);
            size_t index_orig = c_slice.index(j1-cs[1].begin(),j2-cs[2].begin(),i3,i4);
            aux[index_dest] = operator[](index_orig);
            aux.index[index_orig] = index_dest;
          }
        }
      }
    }
  } else if(mask == (MASK_BLOCK3|MASK_BLOCK4)){
    Basis basis(_b3,_b4);
    basis.reorder();

    Basis aux_basis(1);
    aux_basis.reorder();

    aux.resize(_b1,_b2,aux_basis.subspaces(),basis.subspaces());

    for(siter = aux.subspace_begin(); siter != aux.subspace_end(); siter++){
      StateSpace ss = *siter;

      state_slice<T> aux_slice = aux(ss);
      for(int i4 = 0; i4 < ss[4].dim(); i4++){
        State sj = basis(i4+ss[4].begin());
        int j3 = sj.i1;
        int j4 = sj.i2;
        StateSpace cs = get_qn_space(ss[1].qn(),ss[2].qn(),sj.qn1(),sj.qn2());
        cstate_slice<T> c_slice = operator()(cs);

        for(int i1 = 0; i1 < ss[1].dim(); i1++){
          for(int i2 = 0; i2 < ss[2].dim(); i2++){
            size_t index_dest = aux_slice.index(i1,i2,0,i4);
            size_t index_orig = c_slice.index(i1,i2,j3-cs[3].begin(),j4-cs[4].begin());
            aux[index_dest] = operator[](index_orig);
            aux.index[index_orig] = index_dest;
          }
        }
      }
    }
  } else if(mask == (MASK_BLOCK1|MASK_BLOCK3)){
    Basis basis(_b1,_b3);
    basis.reorder();

    Basis aux_basis(1);
    aux_basis.reorder();

    aux.resize(basis.subspaces(),_b2,aux_basis.subspaces(),_b4);

    for(siter = aux.subspace_begin(); siter != aux.subspace_end(); siter++){
      StateSpace ss = *siter;

      state_slice<T> aux_slice = aux(ss);
      for(int i1 = 0; i1 < ss[1].dim(); i1++){
        State sj = basis(i1+ss[1].begin());
        int j1 = sj.i1;
        int j3 = sj.i2;         StateSpace cs = get_qn_space(sj.qn1(),ss[2].qn(),sj.qn2(),ss[4].qn());
        cstate_slice<T> c_slice = operator()(cs);

        for(int i2 = 0; i2 < ss[2].dim(); i2++){
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            size_t index_dest = aux_slice.index(i1,i2,0,i4);
            size_t index_orig = c_slice.index(j1-cs[1].begin(),i2,j3-cs[3].begin(),i4);
            aux[index_dest] = operator[](index_orig);
            aux.index[index_orig] = index_dest;
          }
        }
      }
    }
  } else if(mask == (MASK_BLOCK1|MASK_BLOCK4)){
    Basis basis(_b1,_b4);
    basis.reorder();

    Basis aux_basis(1);
    aux_basis.reorder();

    aux.resize(basis.subspaces(),_b2,_b3,aux_basis.subspaces());

    for(siter = aux.subspace_begin(); siter != aux.subspace_end(); siter++){
      StateSpace ss = *siter;

      state_slice<T> aux_slice = aux(ss);
      for(int i1 = 0; i1 < ss[1].dim(); i1++){
        State sj = basis(i1+ss[1].begin());
        int j1 = sj.i1;
        int j4 = sj.i2;
        StateSpace cs = get_qn_space(sj.qn1(),ss[2].qn(),ss[3].qn(),sj.qn2());
        cstate_slice<T> c_slice = operator()(cs);

        for(int i2 = 0; i2 < ss[2].dim(); i2++){
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            size_t index_dest = aux_slice.index(i1,i2,i3,0);
            size_t index_orig = c_slice.index(j1-cs[1].begin(),i2,i3,j4-cs[4].begin());
            aux[index_dest] = operator[](index_orig);
            aux.index[index_orig] = index_dest;
          }
        }
      }
    }
  } else if(mask == (MASK_BLOCK2|MASK_BLOCK3)){
    Basis basis(_b2,_b3);
    basis.reorder();

    Basis aux_basis(1);
    aux_basis.reorder();

    aux.resize(_b1,basis.subspaces(),aux_basis.subspaces(),_b4);

    for(siter = aux.subspace_begin(); siter != aux.subspace_end(); siter++){
      StateSpace ss = *siter;

      state_slice<T> aux_slice = aux(ss);
      for(int i2 = 0; i2 < ss[2].dim(); i2++){
        State sj = basis(i2+ss[2].begin());
        int j2 = sj.i1;
        int j3 = sj.i2;
        StateSpace cs = get_qn_space(ss[1].qn(),sj.qn1(),sj.qn2(),ss[4].qn());
        cstate_slice<T> c_slice = operator()(cs);

        for(int i1 = 0; i1 < ss[1].dim(); i1++){
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            size_t index_dest = aux_slice.index(i1,i2,0,i4);
            size_t index_orig = c_slice.index(i1,j2-cs[2].begin(),j3-cs[3].begin(),i4);
            aux[index_dest] = operator[](index_orig);
            aux.index[index_orig] = index_dest;
          }
        }
      }
    }


  }
  return aux;
}

////////////////////////////////////////////////////////////////////
// decondense:
// this function rebuilds a state vector from a condensed vector.
////////////////////////////////////////////////////////////////////
template<class T>
VectorState<T>
VectorState<T>::decondense (int mask, const VectorState<T> &orig) const
{
  VectorState<T> aux(orig);

  for(int i = 0; i < Vector<T>::size(); i++)
    aux[i] = operator[](index[i]);

/*
  typename VectorState<T>::iterator siter;

  if(mask == (MASK_BLOCK1|MASK_BLOCK2)){
    Basis basis(aux.b1(),aux.b2());
    basis.reorder();

    for(siter = aux.subspace_begin(); siter != aux.subspace_end(); siter++){
      StateSpace ss(*siter);

      state_slice<T> aux_slice = aux(ss);

      QN qn12 = ss[1].qn() + ss[2].qn();
      StateSpace cs = get_qn_space(qn12,QN(),ss[3].qn(),ss[4].qn());
      cstate_slice<T> s_slice = operator()(cs);

      for(int i1 = 0; i1 < cs[1].dim(); i1++){

        if(basis(i1+cs[1].begin()).qn1() != ss[1].qn() || basis(i1+cs[1].begin()).qn2() != ss[2].qn())
          continue;

        size_t j1 = basis(i1+cs[1].begin()).i1 - ss[1].begin();
        size_t j2 = basis(i1+cs[1].begin()).i2 - ss[2].begin();

        for(int i3 = 0; i3 < ss[3].dim(); i3++){
          for(int i4 = 0; i4 < ss[4].dim(); i4++){

            aux_slice(j1,j2,i3,i4) = s_slice(i1,0,i3,i4);

//            cout << aux_slice.index(j1,j2,i3,i4) << " " << j1+ss[1].begin() << " " << j2+ss[2].begin() << " " << i3 + ss[3].begin() << " " << i4 + ss[4].begin() << "|" << i1+cs[1].begin() << " " << i3+ss[3].begin() << " " << i4 + ss[4].begin() << endl; 
          }

        }
      } 
    }
  } else if(mask == (MASK_BLOCK3|MASK_BLOCK4)){ 
    Basis basis(aux.b3(),aux.b4());
    basis.reorder();

    for(siter = aux.subspace_begin(); siter != aux.subspace_end(); siter++){
      StateSpace ss(*siter);

      state_slice<T> aux_slice = aux(ss);

      QN qn34 = ss[3].qn() + ss[4].qn();
      StateSpace cs = get_qn_space(ss[1].qn(),ss[2].qn(),QN(),qn34);
      cstate_slice<T> s_slice = operator()(cs);

      for(int i4 = 0; i4 < cs[4].dim(); i4++){
        if(basis(i4+cs[4].begin()).qn1() != ss[3].qn() || basis(i4+cs[4].begin()).qn2() != ss[4].qn())
          continue;

        size_t j3 = basis(i4+cs[4].begin()).i1 - ss[3].begin();
        size_t j4 = basis(i4+cs[4].begin()).i2 - ss[4].begin();

        for(int i1 = 0; i1 < ss[1].dim(); i1++){
          for(int i2 = 0; i2 < ss[2].dim(); i2++){

            aux_slice(i1,i2,j3,j4) = s_slice(i1,i2,0,i4);

//            cout << aux_slice.index(i1,i2,j3,j4) << " " << i1+ss[1].begin() << " " << i2+ss[2].begin() << " " << j3 + ss[3].begin() << " " << j4 + ss[4].begin() << "|" << i1+ss[1].begin() << " " << i2+ss[2].begin() << " " << i4 + cs[4].begin() << endl; 
          }

        }
      } 
    }
  } else if(mask == (MASK_BLOCK2|MASK_BLOCK3)){ 
    Basis basis(aux.b2(),aux.b3());
    basis.reorder();

    for(siter = aux.subspace_begin(); siter != aux.subspace_end(); siter++){
      StateSpace ss(*siter);

      state_slice<T> aux_slice = aux(ss);

      QN qn23 = ss[2].qn() + ss[3].qn();
      StateSpace cs = get_qn_space(ss[1].qn(),qn23,QN(),ss[4].qn());
      cstate_slice<T> s_slice = operator()(cs);

      for(int i2 = 0; i2 < cs[2].dim(); i2++){
        if(basis(i2+cs[2].begin()).qn1() != ss[2].qn() || basis(i2+cs[2].begin()).qn2() != ss[3].qn())
          continue;

        size_t j2 = basis(i2+cs[2].begin()).i1 - ss[2].begin();
        size_t j3 = basis(i2+cs[2].begin()).i2 - ss[3].begin();

        for(int i1 = 0; i1 < ss[1].dim(); i1++){
          for(int i4 = 0; i4 < ss[4].dim(); i4++){

            aux_slice(i1,j2,j3,i4) = s_slice(i1,i2,0,i4);

//            cout << aux_slice.index(i1,i2,j3,j4) << " " << i1+ss[1].begin() << " " << i2+ss[2].begin() << " " << j3 + ss[3].begin() << " " << j4 + ss[4].begin() << "|" << i1+ss[1].begin() << " " << i2+ss[2].begin() << " " << i4 + cs[4].begin() << endl; 
          }

        }
      } 
    }


  }
*/

  return aux;
}

////////////////////////////////////////////////////////////////////////
// new_seed:
// rotate ground state of the previous iteration to build seed
// for the new lanczos run.
////////////////////////////////////////////////////////////////////////
template<class T>
void new_seed(const VectorState<T>& v, VectorState<T>& res,
            const BMatrix<T>& rho1, const BMatrix<T> &rho2,
            const Basis& basis1, const Basis& basis2, int position)
{
  VectorState<T> _v = v;

  typename VectorState<T>::const_iterator siter;

  if(position == LEFT){

    Basis aux_basis(1);
    aux_basis.reorder();
    VectorState<T> aux(res.b1(),aux_basis.subspaces(),_v.b3(),_v.b4(),v.qn(),v.qn_mask());
    res = T(0);

    for(siter = aux.subspace_begin(); siter != aux.subspace_end(); siter++){
      StateSpace ss = (*siter);

      const SubMatrix<T> *_block = rho1.block(ss[1].qn());
      if(!_block) continue;
      const SubMatrix<T> &block = *_block;

      Range col_range = block.col_range();
      Range row_range = block.row_range();

      state_slice<T> aux_slice = aux(ss);

      for(int row = row_range.begin(); row <= row_range.end(); row++){
        State sj = basis1(row);
        int j1 = sj.i1;
        int j2 = sj.i2;
        size_t irow = row - row_range.begin();

        StateSpace cs = _v.get_qn_space(sj.qn1(),sj.qn2(),ss[3].qn(),ss[4].qn());
        state_slice<T> s_slice = _v(cs);

        for(int i3 = 0; i3 < ss[3].dim(); i3++){
          for(int i4 = 0; i4 < ss[4].dim(); i4++){

            slice_iter<T> sub_aux = aux_slice(slice(0,aux_slice.size1(),1),0,i3,i4);
            T x = s_slice(j1-cs[1].begin(),j2-cs[2].begin(),i3,i4);

            for(int col = col_range.begin(); col <= col_range.end(); col++){
              size_t icol = col - col_range.begin();
 
              sub_aux(icol) += conj(block(icol,irow))*x;
            }
          }

        }
      }

    }

    for(siter = res.subspace_begin(); siter != res.subspace_end(); siter++){
      StateSpace ss(*siter);

      QN qn34 = ss[3].qn() + ss[4].qn();
      const SubMatrix<T> *_block = rho2.block(qn34);
      if(!_block) continue;
      const SubMatrix<T> &block = *_block;

      Range col_range = block.col_range();
      Range row_range = block.row_range();

      StateSpace cs = aux.get_qn_space(ss[1].qn(),QN(),ss[2].qn(),qn34);
      state_slice<T> aux_slice = aux(cs);
      state_slice<T> res_slice = res(ss);

      for(int row = row_range.begin(); row <= row_range.end(); row++){
        size_t j3 = basis2(row).i1;
        size_t j4 = basis2(row).i2;
        size_t irow = row - row_range.begin();

        if(basis2(row).qn1() != ss[3].qn() || basis2(row).qn2() != ss[4].qn())
          continue;

        for(int i1 = 0; i1 < ss[1].dim(); i1++){
          for(int i2 = 0; i2 < ss[2].dim(); i2++){

            slice_iter<T> sub_aux = aux_slice(i1,0,i2,slice(0,aux_slice.size4(),1));
            gslice_iter<T> sub_res = res_slice(i1,i2,slice(0,res_slice.size3(),1),slice(0,res_slice.size4(),1));
 
            T &val = sub_res(j3-ss[3].begin(),j4-ss[4].begin());

            for(int col = col_range.begin(); col <= col_range.end(); col++){
              size_t icol = col - col_range.begin();

              val += block(icol,irow)*sub_aux(icol);
            }
          }

        }
      }
 
    }
  } else { // position == RIGHT

    Basis aux_basis(1);
    aux_basis.reorder();
    VectorState<T> aux(_v.b1(),_v.b2(),aux_basis.subspaces(),res.b4(),v.qn(),v.qn_mask());
    res = T(0);

    for(siter = aux.subspace_begin(); siter != aux.subspace_end(); siter++){
      StateSpace ss(*siter);

      const SubMatrix<T> *_block = rho1.block(ss[4].qn());
      if(!_block) continue;
      const SubMatrix<T> &block = *_block;

      Range col_range = block.col_range();
      Range row_range = block.row_range();

      state_slice<T> aux_slice = aux(ss);

      for(int row = row_range.begin(); row <= row_range.end(); row++){
        State sj = basis1(row);
        int j3 = sj.i1;
        int j4 = sj.i2;
        size_t irow = row - row_range.begin();

        StateSpace cs = _v.get_qn_space(ss[1].qn(),ss[2].qn(),sj.qn1(),sj.qn2());
        state_slice<T> s_slice = _v(cs);

        for(int i1 = 0; i1 < ss[1].dim(); i1++){
          for(int i2 = 0; i2 < ss[2].dim(); i2++){

            slice_iter<T> sub_aux = aux_slice(i1,i2,0,slice(0,aux_slice.size4(),1));
            T x = s_slice(i1,i2,j3-cs[3].begin(),j4-cs[4].begin());

            for(int col = col_range.begin(); col <= col_range.end(); col++){
              size_t icol = col - col_range.begin();

              sub_aux(icol) += conj(block(icol,irow))*x;
            }
          }

        }
      }

    }

    for(siter = res.subspace_begin(); siter != res.subspace_end(); siter++){
      StateSpace ss(*siter);

      QN qn12 = ss[1].qn() + ss[2].qn();
      const SubMatrix<T> *_block = rho2.block(qn12);
      if(!_block) continue;
      const SubMatrix<T> &block = *_block;

      Range col_range = block.col_range();
      Range row_range = block.row_range();

      StateSpace cs = aux.get_qn_space(qn12,ss[3].qn(),QN(),ss[4].qn());
      state_slice<T> aux_slice = aux(cs);
      state_slice<T> res_slice = res(ss);

      for(int row = row_range.begin(); row <= row_range.end(); row++){
        size_t j1 = basis2(row).i1;
        size_t j2 = basis2(row).i2;
        size_t irow = row - row_range.begin();

        if(basis2(row).qn1() != ss[1].qn() || basis2(row).qn2() != ss[2].qn())
          continue;

        for(int i3 = 0; i3 < ss[3].dim(); i3++){
          for(int i4 = 0; i4 < ss[4].dim(); i4++){

            slice_iter<T> sub_aux = aux_slice(slice(0,aux_slice.size1(),1),i3,0,i4);
            gslice_iter<T> sub_res = res_slice(slice(0,res_slice.size1(),1),slice(0,res_slice.size2(),1),i3,i4);

            T &val = sub_res(j1-ss[1].begin(),j2-ss[2].begin());

            for(int col = col_range.begin(); col <= col_range.end(); col++){
              size_t icol = col - col_range.begin();

              val += block(icol,irow)*sub_aux(icol);
            }
          }

        }
      }
 
    }

  }

  T m = sqrt(real(product(res,res)));
#ifdef DMTK_DEBUG
  cout << "NORM = " << m << endl;
#endif // DMTK_DEBUG
//  res /= m;
 
}

////////////////////////////////////////////////////////

} // namespace dmtk

#endif // __DMTK_STATE_H__
