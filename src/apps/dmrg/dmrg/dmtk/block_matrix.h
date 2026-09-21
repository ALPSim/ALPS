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

#ifndef __DMTK_BLOCK_MATRIX_H__
#define __DMTK_BLOCK_MATRIX_H__

#include <iosfwd>
#include <vector>
#include <complex>
#include "conj.h"
#include "matrix.h"
#include "qn.h"
#include "subspace.h"
#include "basis.h"
#ifdef WITH_LAPACK
#include "lapack_interface.h"
#endif
#include "lanczos.cc"
#include "util.h" 

namespace dmtk
{

template<class T>
class SubMatrix: public Matrix<T>
{
  private:
    QN _qn;
    SubSpace _row_range;
    SubSpace _col_range;
  public:
    SubMatrix(): Matrix<T>(), _qn(0), _row_range(QN(),-1,-1), _col_range(QN(),-1,-1) {}
    SubMatrix(QN qn, SubSpace col_range, SubSpace row_range): 
      Matrix<T>(col_range.size(), row_range.size()), _qn(qn), _row_range(row_range), _col_range(col_range) {}
    SubMatrix(const SubMatrix<T> &s): 
      Matrix<T>(s), _qn(s._qn), _row_range(s._row_range), _col_range(s._col_range) {}

//  Asignment

    SubMatrix& operator=(const SubMatrix<T> &sm)
      { Matrix<T>::operator=(sm); _qn = sm._qn; _col_range = sm._col_range; _row_range = sm._row_range; return *this; }
    SubMatrix& operator=(const Matrix<T> &m)
      { Matrix<T>::operator=(m); return *this; }
    SubMatrix& operator=(const Sparse<T>& s)
      { Matrix<T>::operator=(s); return *this; }

    SubMatrix& operator=(cgslice_iter<T> s)
      { Matrix<T>::operator=(s); return *this; }
    SubMatrix& operator=(gslice_iter<T> s)
      { Matrix<T>::operator=(s); return *this; }
    SubMatrix& operator=(const T *v)
      { Matrix<T>::operator=(v); return *this; }
    SubMatrix& operator=(const T& v)
      { Matrix<T>::operator=(v); return *this; }
    template<class Expr>
    SubMatrix& operator=(const IterExpr<T,Expr>&expr)
      { Matrix<T>::operator=(expr); return *this; }

    SubMatrix& operator+=(const Matrix<T> &m)
      { Matrix<T>::operator+=(m); return *this; }
    SubMatrix& operator-=(const Matrix<T> &m)
      { Matrix<T>::operator-=(m); return *this; }
    SubMatrix& operator*=(const Matrix<T> &m)
      { Matrix<T>::operator*=(m); return *this; }
    SubMatrix& operator/=(const Matrix<T> &m)
      { Matrix<T>::operator/=(m); return *this; }

    SubMatrix& operator+=(const T &v)
      { Matrix<T>::operator+=(v); return *this; }
    SubMatrix& operator-=(const T &v)
      { Matrix<T>::operator-=(v); return *this; }
    SubMatrix& operator*=(const T &v)
      { Matrix<T>::operator*=(v); return *this; }
    SubMatrix& operator/=(const T &v)
      { Matrix<T>::operator/=(v); return *this; }

    template<class Expr>
    SubMatrix& operator+=(const IterExpr<T,Expr>&expr)
      { Matrix<T>::operator+=(expr); return *this; }
    template<class Expr>
    SubMatrix& operator-=(const IterExpr<T,Expr>&expr)
      { Matrix<T>::operator-=(expr); return *this; }
    template<class Expr>
    SubMatrix& operator*=(const IterExpr<T,Expr>&expr)
      { Matrix<T>::operator*=(expr); return *this; }
    template<class Expr>
    SubMatrix& operator/=(const IterExpr<T,Expr>&expr)
      { Matrix<T>::operator/=(expr); return *this; }

//  Compare

    bool operator==(const SubMatrix<T> &m)
      {
        if(_qn == m._qn && _row_range == m._row_range && _col_range == m._row_range) return true;
        return false;
      }

//  Methods

    SubMatrix<T>& resize(const Range &col_range, const Range &row_range)
      {
        Matrix<T>::resize(col_range.size(), row_range.size());
        _col_range = col_range;
        _row_range = row_range;
        return *this;
      }
    SubMatrix<T>& resize(size_t cols, size_t rows)
      {
        Matrix<T>::resize(cols, rows);
        _col_range = Range(0,cols-1);
        _row_range = Range(0,rows-1);
        return *this;
      }
    SubMatrix<T>& reshape(const Range &col_range, const Range &row_range)
      {
        Matrix<T>::reshape(col_range.size(), row_range.size());
        _col_range = col_range;
        _row_range = row_range;
        return *this;
      }
    SubMatrix<T>& reshape(size_t cols, size_t rows)
      {
        Matrix<T>::reshape(cols, rows);
        _col_range = Range(0,cols-1);
        _row_range = Range(0,rows-1);
        return *this;
      }

    SubSpace col_range() const { return _col_range; }
    SubSpace row_range() const { return _row_range; }
    QN qn() const { return _qn; }
    QN& qn() { return _qn; }

//  Streams

    void read(std::istream& s)
    {
      _qn.read(s);
      _row_range.read(s); 
      _col_range.read(s); 
      Matrix<T>::read(s); 
    }

    void write(std::ostream& s) const
    {
      _qn.write(s);
      _row_range.write(s); 
      _col_range.write(s); 
      Matrix<T>::write(s); 
    }

};

template <class T>
class BMatrix: public std::vector<SubMatrix<T> >
{
  private:
    PackedBasis _subspace;
    Vector<int> _hash_index;
    Vector<int> _dqn;
    Vector<int> _qn_min;

#ifdef USE_HASH
    int _get_index(const QN& qn) const
      {
        int index = 0;
        for(int i = 0, j = 0; i < QN::QN_LAST; i++) {
          if(QN::get_qn_mask() & (1 << i)){
            index += (qn[i]-_qn_min[j])*_dqn[j];
            j++;
          }
        }
        return index; 
      }

    void init_hash()
      {
        _dqn.resize(QN::QN_LAST);
        _qn_min.resize(QN::QN_LAST);

        PackedBasis::iterator iter;
        iter = _subspace.begin();
        Vector<int> min(QN::QN_LAST,999999);
        iter = _subspace.end();
        iter--;
        Vector<int> max(QN::QN_LAST,-999999);
        int jmax;
        for(iter = _subspace.begin(); iter != _subspace.end(); iter++){
          QN &qn = (*iter).qn();
          for(int i = 0, j = 0; i < QN::QN_LAST; i++){
            if(QN::get_qn_mask() & (1 << i)){
              int qni = qn[i];
              if(qni < min[j]) min[j] = qni;
              if(qni > max[j]) max[j] = qni;
              j++;
            }
            jmax = j;
          }
        }
        int dim = 1;
        Vector<int> dqn(QN::QN_LAST);
        for(int j = 0; j < jmax; j++) {
          dqn[j] = max[j]-min[j]+1;
          _qn_min[j] = min[j];
          dim *= dqn[j];
        }
        _hash_index.resize(dim);
        
        dim = 1;
        _dqn = 1;
        for(int j = jmax-2; j >= 0; j--) {
          dim *= dqn[j+1];
          _dqn[j] = dim;
        }

        iterator viter;
        _hash_index = -1;
        int idx = 0;
        for(viter = _V::begin(); viter != _V::end(); viter++){
          _hash_index[_get_index((*viter).qn())] = idx++;
        }
      }
#else // !USE_HASH
    int _get_index(const QN& qn) const { return -1; }
    void init_hash() {};
#endif // USE_HASH

    BMatrix& init(Basis &b)
      {
        std::vector<SubMatrix<T> >::clear();
        b.reorder();
        _subspace = b.subspaces();
        init_hash();
        return *this;
      }

  public:
    typedef typename std::vector<SubMatrix<T> > _V;
    typedef typename std::vector<SubMatrix<T> >::iterator iterator;
    typedef typename std::vector<SubMatrix<T> >::const_iterator const_iterator;

    BMatrix(): _V() { _dqn.resize(QN::QN_LAST); _qn_min.resize(QN::QN_LAST); }
    BMatrix(const BMatrix<T>& m):_V(m),_subspace(m._subspace), _hash_index(m._hash_index), _dqn(m._dqn), _qn_min(m._qn_min) {}
    BMatrix(const PackedBasis& b): _subspace(b) { init_hash(); }
    BMatrix(const Basis& b) { Basis _b(b); init(_b); }
    BMatrix(const Basis& b1, const Basis &b2) { Basis b(b1,b2); init(b); }

    BMatrix& repack(const Basis &basis)
      { _subspace = basis.subspaces(); init_hash(); return *this; }
    BMatrix& repack(const PackedBasis &basis)
      { _subspace = basis; init_hash(); return *this; }

#ifdef USE_HASH
    BMatrix& push_back(const SubMatrix<T> &o)
      {
        _V::push_back(o);
        _hash_index[_get_index(o.qn())] = _V::size()-1;
        return *this;
      }

    iterator erase(iterator iter)
      {
        QN qn = (*iter).qn();
        iterator biter;

        _hash_index[_get_index(qn)] = -1;
        biter = iter;
        biter++;
        for(; biter != _V::end(); biter++){
          _hash_index[_get_index((*biter).qn())]--;
        }

        iterator res_iter = _V::erase(iter);
        return res_iter;
      } 
#endif  // USE_HASH

    BMatrix& operator=(const BMatrix& m)
      {
        _V::operator=(m);

        _subspace = m._subspace;
        _hash_index = m._hash_index;
        _dqn = m._dqn;
        _qn_min = m._qn_min;

        return *this;
      }

/* ORIGINAL
    BMatrix& operator+=(const BMatrix& v)
      {
        iterator iter;
        const_iterator citer;
        for(iter = _V::begin(), citer = v.begin(); iter != _V::end() && citer != v.end(); iter++, citer++)
          {
            SubMatrix<T> &m = (*iter);
            m += (*citer);
          }
        return *this;
      }
*/
    BMatrix& operator+=(const BMatrix& v)
      {
        iterator iter;
        const_iterator citer;
        for(iter = _V::begin(); iter != _V::end(); iter++)
          {
          for(citer = v.begin(); citer != v.end(); citer++)
            {
              if((*citer).qn() == (*iter).qn()){
                SubMatrix<T> &m = (*iter);
                m += (*citer);
              }
            }
          }
        return *this;
      }


    BMatrix operator*(const T& v) const
      {
        BMatrix<T> aux(*this);
        iterator iter;
        for(iter = aux.begin(); iter != aux.end(); iter++)
          {
            SubMatrix<T> &m = (*iter);
            m *= v;
          }

        return aux;
      }

    BMatrix& operator=(const T& v)
      {
        iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++)
          {
            SubMatrix<T> &m = (*iter);
            m = v;
          }
        return *this;
      }

    const SubMatrix<T>* block(const QN &qn) const
      {
        const SubMatrix<T> *sm = get_block(qn);
        return sm;
      }     

    SubMatrix<T>* block(const QN &qn)
      {
        SubMatrix<T> *sm = get_block(qn);
        return sm;
      }     

    SubSpace subspace(const QN &qn) const
      {
        return(_subspace(qn));
      }     

    SubSpace subspace(size_t index) const
      {
        return _subspace[index];        
      }     

    const PackedBasis& subspaces() const { return _subspace; }

    SubMatrix<T> *operator[](size_t index)
      {
        SubMatrix<T> &sm = _V::operator[](index);
        return &sm;
      }

    const SubMatrix<T> *operator[](size_t index) const
      {
        const SubMatrix<T> &sm = _V::operator[](index);
        return &sm;
      }

    T operator() (size_t col, size_t row) const
      {
        const_iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++){
          const SubMatrix<T> &b = (*iter);
          if(col >= b.col_range().begin() && col <= b.col_range().end() &&
             row >= b.row_range().begin() && row <= b.row_range().end())
               return b(col-b.col_range().begin(),row-b.row_range().begin());
        }
//        return T(-99999);
        return T(0);
      }

    PackedBasis::iterator subspace_begin() { return _subspace.begin(); }
    PackedBasis::const_iterator subspace_begin() const { return _subspace.begin(); }
    PackedBasis::iterator subspace_end() { return _subspace.end(); }
    PackedBasis::const_iterator subspace_end() const { return _subspace.end(); }

    SubMatrix<T> *get_block(const QN& qn);
    const SubMatrix<T> *get_block(const QN& qn) const;
    int get_block_index(const QN& qn) const;

    BMatrix& diagonalize(Vector<double>& ev);
    BMatrix& diagonalize(Vector<float>& ev);

    size_t dim() const
     {
        PackedBasis::const_iterator iter;
        size_t d = 0;
        for(iter = _subspace.begin(); iter != _subspace.end(); iter++)
          d += (*iter).dim();
        return d;        
     }

    BMatrix& resize(const Basis& b) 
      { if(_V::size() != 0) _V::clear(); Basis _b(b); init(_b); return *this; }
    BMatrix& resize(const Basis& b1, const Basis &b2) 
      { if(_V::size() != 0) _V::clear(); Basis b(b1,b2); init(b); return *this; }

//  Streams

    void read(std::istream& s)
    {
      size_t l;
      s.read((char *)&l, sizeof(size_t));
      _V::resize(l);

      s.read((char *)&l, sizeof(size_t));
      _subspace.resize(l);

      PackedBasis::iterator siter;
      for(siter = _subspace.begin(); siter != _subspace.end(); siter++){
        (*siter).read(s);
//        cout << (*siter).qn().n() << " " << (*siter).qn().sz() << endl;
      }

      _hash_index.read(s);
      _dqn.read(s);
      _qn_min.read(s);

      iterator iter;
      for(iter = _V::begin(); iter != _V::end(); iter++){
        SubMatrix<T> &sm = (*iter);
        sm.read(s);
      }
    }

    void write(std::ostream& s) const
    {
      size_t l = _V::size();
      s.write((const char *)&l, sizeof(size_t));

      l = _subspace.size();
      s.write((const char *)&l, sizeof(size_t));

      PackedBasis::const_iterator siter;
      for(siter = _subspace.begin(); siter != _subspace.end(); siter++)
        (*siter).write(s);

      _hash_index.write(s);
      _dqn.write(s);
      _qn_min.write(s);

      const_iterator iter;
      for(iter = _V::begin(); iter != _V::end(); iter++){
        const SubMatrix<T> &sub = (*iter);
        sub.write(s);
      }
    }

};

//////////////////////////////////////////////////////////////////////
template<class T>
SubMatrix<T>*
BMatrix<T>::get_block(const QN &qn)
{
#ifdef USE_HASH
  int qn_index = _get_index(qn);
  if(qn_index < 0 || qn_index >= _hash_index.size()) return NULL;
  int idx = _hash_index[qn_index];
  if(idx >= 0 && idx < size()) 
    return operator[](idx);
  return NULL;
#endif // USE_HASH
/*
  for(int i = 0; i < this->size(); i++){
    SubMatrix<T> *_s = operator[](i);
    cout  << "HOLA BLOCK " << _s->qn() << endl;
//    if(_s->qn() == qn) return _s; 
  }
*/

  int origin = 0;
  int end = _V::size() - 1;
  int index = -1;

  while(origin <= end){
    int index_old = index;
    index = (origin + end) / 2;

    SubMatrix<T> *_s = operator[](index);

    if(_s->qn() == qn) return _s; 
    if(_s->qn() > qn)
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

/*
cout << "HOLA BLOCK NOT FOUND " << qn[0] << endl;
  for(int i = 0; i < this->size(); i++){
    const SubMatrix<T> *_s = operator[](i);
    if(_s->qn() == qn) cout << qn[0] << endl; 
  }
*/
  return NULL;
}

template<class T>
const SubMatrix<T>*
BMatrix<T>::get_block(const QN &qn) const
{
#ifdef USE_HASH
  int qn_index = _get_index(qn);
  if(qn_index < 0 || qn_index >= _hash_index.size()) return NULL;
  int idx = _hash_index[qn_index];
  if(idx >= 0 && idx < size()) 
    return operator[](idx);
  return NULL;
#endif // USE_HASH
/*
  for(int i = 0; i < this->size(); i++){
    const SubMatrix<T> *_s = operator[](i);
    if(_s->qn() == qn) return _s; 
  }
*/

  int origin = 0;
  int end = _V::size() - 1;
  int index = -1;

  while(origin <= end){
    int index_old = index;
    index = (origin + end) / 2;

    const SubMatrix<T> *_s = operator[](index);

    if(_s->qn() == qn) return _s; 
    if(_s->qn() > qn)
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

/*
cout << "HOLA BLOCK NOT FOUND " << qn[0] << endl;
  for(int i = 0; i < this->size(); i++){
    const SubMatrix<T> *_s = operator[](i);
    if(_s->qn() == qn) cout << qn[0] << endl; 
  }
*/
  return NULL;
}

template<class T>
int
BMatrix<T>::get_block_index(const QN &qn) const
{
#ifdef USE_HASH
  int qn_index = _get_index(qn);
  if(qn_index < 0 || qn_index >= _hash_index.size()) return NULL;
  int idx = _hash_index[qn_index];
  if(idx >= 0 && idx < size()) 
    return idx;
  return -1;
#endif // USE_HASH
/*
  for(int i = 0; i < this->size(); i++){
    const SubMatrix<T> *_s = operator[](i);
    if(_s->qn() == qn) return i; 
  }
*/

  int origin = 0;
  int end = _V::size() - 1;
  int index = -1;

  while(origin <= end){
    int index_old = index;
    index = (origin + end) / 2;

    const SubMatrix<T> *_s = operator[](index);

    if(_s->qn() == qn) return index; 
    if(_s->qn() > qn)
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

  return -1;
}

//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////

template<>
BMatrix<float>& 
BMatrix<float>::diagonalize(Vector<float>& _ev)
{
  Vector<float> ev(this->dim()); // eigenvalues

  iterator iter;
  for(iter = begin(); iter != end(); iter++){

#ifdef WITH_LAPACK
    Matrix<float> &b = (*iter);

    Vector<float> d(b.rows());
    Vector<float> work(3*b.rows());
    int info = 0;

    FORTRAN_ID(ssyev)('V','U',b.rows(),b.array(),b.rows(),
           d.array(),work.array(),work.size(),info);

    if(info != 0) cout << "*** ERROR in dsyev: info != 0 (failed to converge)\n";

    ev((*iter).row_range()) = d(Range(0,b.rows()-1));

#endif // WITH LAPACK

  }

  _ev.resize(ev.size());
  for(int i = 0; i < ev.size(); i++) _ev[i] = ev[i];
  return *this;
}


template<>
BMatrix<double>& 
BMatrix<double>::diagonalize(Vector<double>& ev)
{
  ev.resize(this->dim()); // eigenvalues

  iterator iter;
  for(iter = begin(); iter != end(); iter++){

#ifdef WITH_LAPACK
    Matrix<double> &b = (*iter);

    Vector<double> d(b.rows());
    Vector<double> work(3*b.rows());
    int info;

    FORTRAN_ID(dsyev)('V','U',b.rows(),b.array(),b.rows(),
           d.array(),work.array(),work.size(),info);

    if(info != 0) cout << "*** ERROR in dsyev: info != 0 (failed to converge)\n";

    ev((*iter).row_range()) = d(Range(0,b.rows()-1));

#endif // WITH LAPACK

  }

  return *this;
}

template<>
BMatrix<std::complex<double> >& 
BMatrix<std::complex<double> >::diagonalize(Vector<double>& ev)
{
  ev.resize(this->dim()); // eigenvalues

  iterator iter;
  for(iter = begin(); iter != end(); iter++){

    Matrix<std::complex<double> > &b = (*iter);

#ifdef WITH_LAPACK

    Vector<double> d(b.rows());
    Vector<std::complex<double> > work(2*b.rows());
    Vector<double> rwork(3*b.rows()-2);

    int info;

    FORTRAN_ID(zheev)('V','U',b.rows(),b.array(),b.rows(),
           d.array(),work.array(),work.size(),rwork.array(),info);

    if(info != 0) cout << "*** ERROR in zheev: info != 0 (failed to converge)\n";

    ev((*iter).row_range()) = d(Range(0,b.rows()-1));

#endif // WITH_LAPACK

  }
  return *this;
}


//////////////////////////////////////////////////////////////////////
template<class T>
BMatrix<T>&
cexp(const BMatrix<T> &bm)
{
  BMatrix<T> aux(bm);
  typename BMatrix<T>::const_iterator citer;
  typename BMatrix<T>::iterator iter;
  for(citer = bm.begin(), iter = aux.begin(); citer != bm.end(); citer++, iter++){
    const Matrix<T> & m = (*citer);
    Matrix<T> &xm = (*iter);
    xm = cexp(m);
  }
  return aux;
} 

template<class T>
BMatrix<T>&
exp(const BMatrix<T> &bm)
{
  BMatrix<T> aux(bm);
  typename BMatrix<T>::const_iterator citer;
  typename BMatrix<T>::iterator iter;
  for(citer = bm.begin(), iter = aux.begin(); citer != bm.end(); citer++, iter++){
    const Matrix<T> & m = (*citer);
    Matrix<T> &xm = (*iter);
    xm = exp(m);
  }
  return aux;
}

template<class T>
void
biorthonormalize(BMatrix<T> &a, BMatrix<T> &b)
{
  BMatrix<T> ab, u, v;
  ab = a;
  ab.clear();
  Vector<double> ev;

  int dim = 0;
  for(int i = 0; i < a.size(); i++){
    SubMatrix<T> &_a = *a[i];
    SubMatrix<T> &_b = *b[i];
    SubMatrix<T> _ab(_a);
    _ab = product(_a.ct(),_b);
    ab.push_back(_ab);
    dim += _ab.rows();
  }

  ev.resize(dim);
  ab.svd(ev, u, v);

  int idx = 0; 
  for(int i = 0; i < a.size(); i++){ 
    Matrix<T> &_a = *a[i];
    Matrix<T> &_b = *b[i];
    Matrix<T> &_u = *u[i];
    Matrix<T> &_v = *v[i];
    for(int j = 0; j < _a.rows(); j++){
      _u.column(j) = _u.column(j)*(1./sqrt(ev(idx)));
      _v.column(j) = _v.column(j)*(1./sqrt(ev(idx)));
      idx++;
    }
    Matrix<T> aux(_a.cols(), _a.rows()); 
    aux = product(_a,_u), _a = aux;
    aux = product(_b,_v), _b = aux;
//    for(int j = 0; j < _a.rows(); j++){
//      cout << "BI-ORTHO " << j << " " << product(_a.column(j),_b.column(j)) << endl;
//    }
  }
}

//////////////////////////////////////////////////////////////////////

template<class T>
void
product(const SubMatrix<T> &block1, const SubMatrix<T>& block2,
        const Matrix<T>& subv, Matrix<T>& subres, Matrix<T> &aux, 
        T coefa = T(1), T coefb = T(0),
        bool hc = false)
{   
  if(!hc){
    aux.reshape(subv.cols(),block2.rows());
    matrix_matrix_product('N','N',block2,subv,aux);

    matrix_matrix_product('N','T',aux,block1,subres,coefa,coefb);

  } else {
    aux.reshape(subv.cols(),block2.cols());
    matrix_matrix_product('C','N',block2,subv,aux);

    Matrix<T> aux1(block1.cols(),block1.rows());
    aux1 = conj(block1);
    matrix_matrix_product('N','N',aux,aux1,subres,dmtk::conj(coefa),dmtk::conj(coefb));
  }
}

template<class T>
void
product(const SubMatrix<T> &block1, const SubMatrix<T>& block2,
        const Vector<T>& subv, Vector<T>& subres, 
        T coefa = T(1), T coefb = T(0),
        bool hc = false)
{   
  if(!hc){
    Matrix<T> mij(block2.cols(),block1.rows());
    matrix_matrix_product('N','N',block1,block2,mij);
   
    matrix_vector_product('N', mij,subv,subres,coefa,coefb);

  } else {
    Matrix<T> mij(block1.rows(),block2.cols());
    matrix_matrix_product('C','C',block2,block1,mij);
    
    matrix_vector_product('N', mij,subv,subres,dmtk::conj(coefa),coefb);
  }

}

//////////////////////////////////////////////////////////////////////
// Version for condensed blocks
//////////////////////////////////////////////////////////////////////

template<class T>
void
product(const SubMatrix<T> &block1, const SubMatrix<T>& block2,
        cgslice_iter<T>& subv, gslice_iter<T>& subres, Matrix<T> &aux, 
        T coefa = T(1), T coefb = T(0),
        bool hc = false)
{   
  if(!hc){
    aux.reshape(subv.size1(),block2.rows());
    matrix_matrix_product('N','N',block2.rows(),subv.size1(),block2.cols(),block2.array(),block2.rows(),subv.get_pointer(0,0),subv.size2(),aux.array(),aux.rows());

    matrix_matrix_product('N','T',aux.rows(),block1.rows(),aux.cols(),aux.array(),aux.rows(),block1.array(),block1.rows(),&subres(0,0),subres.size2(),coefa,coefb);

  } else {
    aux.reshape(subv.size1(),block2.cols());
    matrix_matrix_product('C','N',block2.cols(),subv.size1(),block2.rows(),block2.array(),block2.rows(),subv.get_pointer(0,0),subv.size2(),aux.array(),aux.rows());

    Matrix<T> aux1(block1.cols(),block1.rows());
    aux1 = conj(block1);
    matrix_matrix_product('N','N',aux.rows(),block1.cols(),aux.cols(),aux.array(),aux.rows(),aux1.array(),block1.rows(),&subres(0,0),subres.size2(),dmtk::conj(coefa),dmtk::conj(coefb));
  }
}

} // namespace dmtk

#endif // __DMTK_BLOCK_MATRIX_H__
