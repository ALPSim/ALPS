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

#ifndef __DMTK_BLOCK_MATRIX_H__
#define __DMTK_BLOCK_MATRIX_H__

#include <iosfwd>
#include <vector>
#include <list>
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
    SubMatrix(): Matrix<T>(), _row_range(QN(),-1,-1), _col_range(QN(),-1,-1) {}
    SubMatrix(QN qn, SubSpace col_range, SubSpace row_range): 
      Matrix<T>(col_range.size(), row_range.size()), _qn(qn), _row_range(row_range), _col_range(col_range) {}

//  Asignment

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
class BMatrix: public std::list<SubMatrix<T> >
{
  private:
    PackedBasis _subspace;
    BMatrix& init(Basis &b)
      {
        std::list<SubMatrix<T> >::clear();
        b.reorder();
        _subspace = b.subspaces();
        return *this;
      }

  public:
    typedef typename std::list<SubMatrix<T> > _V;
    typedef typename std::list<SubMatrix<T> >::iterator iterator;
    typedef typename std::list<SubMatrix<T> >::const_iterator const_iterator;

    BMatrix() {}
    BMatrix(const BMatrix<T>& m):_V(m),_subspace(m._subspace){}
    BMatrix(const Basis& b) { Basis _b(b); init(_b); }
    BMatrix(const PackedBasis& b): _subspace(b) {}
    BMatrix(const Basis& b1, const Basis &b2) { Basis b(b1,b2); init(b); }

    BMatrix& repack(const Basis &basis)
      { _subspace = basis.subspaces(); return *this; }
    BMatrix& repack(const PackedBasis &basis)
      { _subspace = basis; return *this; }

    BMatrix& operator=(const BMatrix& m)
      {
        _V::operator=(m);
        _subspace = m._subspace;
        return *this;
      }
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
// TODO: We could use binary search
        const_iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++){
          if((*iter).qn() == qn) return &(*iter);
        }

        return 0;
      }     

    SubMatrix<T>* block(const QN &qn) 
      {
// TODO: We could use binary search
        iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++){
          if((*iter).qn() == qn) return &(*iter);
        }

        return 0;
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
        typename _V::iterator iter = _V::begin();
        if(_V::size() == 0) return 0; // Sanity check;
        int i;
        for(i = 0; i < index && iter != _V::end(); i++, iter++){}
        if(i != index) return 0;
        return &(*iter);
      }

    const SubMatrix<T> *operator[](size_t index) const
      {
        typename _V::const_iterator iter = _V::begin();
        if(_V::size() == 0) return 0; // Sanity check;
        int i;
        for(i = 0; i < index && iter != _V::end(); i++, iter++){}
        if(i != index) return 0;
        return &(*iter);
      }

    T operator() (size_t col, size_t row) const
      {
        const_iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++){
          const SubMatrix<T> &b(*iter);
          if(col >= b.col_range().begin() && col <= b.col_range().end() &&
             row >= b.row_range().begin() && row <= b.row_range().end())
               return b(col-b.col_range().begin(),row-b.row_range().begin());
        }
        return T(0);
      }

    PackedBasis::iterator subspace_begin() { return _subspace.begin(); }
    PackedBasis::const_iterator subspace_begin() const { return _subspace.begin(); }
    PackedBasis::iterator subspace_end() { return _subspace.end(); }
    PackedBasis::const_iterator subspace_end() const { return _subspace.end(); }

    BMatrix& diagonalize(Vector<double>& ev);

    size_t dim() const
     {
        std::vector<SubSpace>::const_iterator iter;
        size_t d = 0;
        for(iter = _subspace.begin(); iter != _subspace.end(); iter++)
          d += (*iter).dim();
        return d;        
     }

    BMatrix& resize(const Basis& b) 
      { _V::clear(); Basis _b(b); init(_b); return *this; }
    BMatrix& resize(const Basis& b1, const Basis &b2) 
      { _V::clear(); Basis b(b1,b2); init(b); return *this; }

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
      }

      iterator iter;
      for(iter = _V::begin(); iter != _V::end(); iter++){
        SubMatrix<T> &sub = (*iter);
        sub.read(s);
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

      const_iterator iter;
      for(iter = _V::begin(); iter != _V::end(); iter++){
        const SubMatrix<T> &sub = (*iter);
        sub.write(s);
      }
    }

};

//////////////////////////////////////////////////////////////////////

template<>
BMatrix<double>& 
BMatrix<double>::diagonalize(Vector<double>& ev)
{
  ev.resize(this->dim()); // eigenvalues

  iterator iter;
  for(iter = begin(); iter != end(); iter++){

#ifdef WITH_LAPACK

    Matrix<double> &b = *iter;

    Vector<double> d(b.rows());
    Vector<double> work(3*b.rows());
    int info;

    FORTRAN_ID(dsyev)('V','U',b.rows(),b.array(),b.rows(), 
           d.array(),work.array(),work.size(),info);

    if(info != 0) cerr << "*** ERROR in dsyev: info != 0 (failed to converge)\n";

    ev((*iter).row_range()) = d(Range(0,b.rows()-1));

#endif // WITH LAPACK

  }

  return *this;
}

template<>
BMatrix<complex<double> >& 
BMatrix<complex<double> >::diagonalize(Vector<double>& ev)
{
  ev.resize(this->dim()); // eigenvalues

  iterator iter;
  for(iter = begin(); iter != end(); iter++){

    Matrix<complex<double> > &b = *iter;

#ifdef WITH_LAPACK

    Vector<double> d(b.rows());
    Vector<complex<double> > work(2*b.rows());
    Vector<double> rwork(3*b.rows()-2);
    int info;

    FORTRAN_ID(zheev)('V','U',b.rows(),b.array(),b.rows(), 
           d.array(),work.array(),work.size(),rwork.array(),info);

    if(info != 0) cerr << "*** ERROR in zheev: info != 0 (failed to converge)\n";

    ev((*iter).row_range()) = d(Range(0,b.rows()-1));

#endif // WITH_LAPACK

  }
  return *this;
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
    matrix_matrix_product('N','N',aux,aux1,subres,conj(coefa),conj(coefb));

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
    
    matrix_vector_product('N', mij,subv,subres,conj(coefa),coefb);
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
    matrix_matrix_product('C','N',block2.cols(),subv.size1(),block2.cols(),block2.array(),block2.rows(),subv.get_pointer(0,0),subv.size2(),aux.array(),aux.rows());

    Matrix<T> aux1(block1.cols(),block1.rows());
    aux1 = conj(block1);
    matrix_matrix_product('N','N',aux.rows(),block1.cols(),aux.cols(),aux.array(),aux.rows(),aux1.array(),block1.rows(),&subres(0,0),subres.size2(),std::conj(coefa),std::conj(coefb));
  }
}

//////////////////////////////////////////////////////////////
// OLD VERSION
//////////////////////////////////////////////////////////////

/*
template<class T>
void
product(const SubMatrix<T> &block1, const SubMatrix<T>& block2,
        cgslice_iter<T> subv, gslice_iter<T> subres, T coef = T(1),
        bool hc = false)
{
  if(!hc){
    for(int i1 = 0; i1 < block1.rows(); i1++){
      for(int i2 = 0; i2 < block2.rows(); i2++){

        for(int j1 = 0; j1 < block1.cols(); j1++){

          T mij = T(0);
          for(int j2 = 0; j2 < block2.cols(); j2++){
            mij += block2(j2,i2)*subv(j1,j2); 
          }
          subres(i1,i2) += coef*block1(j1,i1)*mij;
        }
      }
    }

  } else {

    for(int i1 = 0; i1 < block1.cols(); i1++){
      for(int i2 = 0; i2 < block2.cols(); i2++){

        for(int j1 = 0; j1 < block1.rows(); j1++){

          T mij = T(0);
          for(int j2 = 0; j2 < block2.rows(); j2++){
            mij += std::conj(block2(i2,j2))*subv(j1,j2);
          }
          subres(i1,i2) += std::conj(coef)*std::conj(block1(i1,j1))*mij;
        }
      }
    }
  }
}
*/
/*
template<class T>
void
product(const SubMatrix<T> &block1, const SubMatrix<T>& block2,
        cslice_iter<T> subv, slice_iter<T> subres, T coef = T(1),
        bool hc = false)
{  
  if(!hc){
    for(int i = 0; i < block1.rows(); i++){
      T v = T(0);
      for(int j = 0; j < block2.cols(); j++){
        T mij = T(0);
        for(int k = 0; k < block2.rows(); k++){
//cout << i << " " << k << " " << block1(k,i) << endl;
//cout << j << " " << k << " " << block2(j,k) << endl;
          mij += block1(k,i)*block2(j,k);
        }
        v += mij * subv(j);
      }
      subres(i) += coef * v;
    }

  } else {
    for(int i = 0; i < block2.cols(); i++){
      T v = T(0);
      for(int j = 0; j < block1.rows(); j++){
        T mij = T(0);
        for(int k = 0; k < block2.rows(); k++){
          mij += std::conj(block2(i,k))*std::conj(block1(k,j));
        }
        v += mij * subv(j);
      }
      subres(i) += std::conj(coef) * v;
    }

  }

}
*/
       
} // namespace dmtk

#endif // __DMTK_BLOCK_MATRIX_H__
