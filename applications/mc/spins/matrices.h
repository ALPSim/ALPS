/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2005 by Matthias Troyer <troyer@comp-phys.org>,
*                       Andreas Streich <astreich@student.ethz.ch>
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

/* $Id$ */

#ifndef ALPS_APPLICATIONS_MC_SPIN_MATRICES_H_
#define ALPS_APPLICATIONS_MC_SPIN_MATRICES_H_

#include <algorithm>
#include <iostream>
#include <alps/parameter.h>

#include "helper.h"

template <class T, int N>
class TinyVector;

template <class T, int N>
class SquareMatrix {
  public:
    /**
      * Constructer which sets all matrix elements to x.
      */
    SquareMatrix(T x=0.) { setMatrix(x); }

    /** 
      * Constructer which parces the given string, evaluates the substrings and 
      * sets the matrix, in row-first manner. The substrings are delimited by 
      * white spaces.
      *
      * \param str the substring to be parced.
      * \param parms contains the values for expressions contained in the 
      *              string.
      */
    SquareMatrix(std::string& str,alps::Parameters parms) 
    { setMatrix(str,parms); }
    
    static const int dim = N;
    static const bool allow_cluster_update = false;
    void setMatrix(std::string&,alps::Parameters); 
    void setMatrix(T x);

    /** 
     * Returns the element at position (i,j) of the matrix. The index bounds
     * are checked, an exception is thrown if the given indexes are not valid. 
     */
    T getElement(int i, int j) const {
      if (((i>=0) && (i < N)) && ((j>=0) && (j<N)))
        return elem_[i*N+j];
      else {
        std::cerr << "Invalid index in TinyMatrix: \n"
                  << "Matrix dimensions are ( " << N << " x " << N << ")\n"
                  << "given indexes are [" << i << "][" << j << "].\n";
        boost::throw_exception(std::out_of_range("invalid index"));
      }
    }
    
    /**
     * Returns the element at position (i,j) of the matrix. The index bounds are
     * NOT checked -- use this function only if you are sure that the indexes 
     * are ok!
     * The version is faster than getElement(int i, int j).
     */
    inline T getElement_nc(int i, int j) const { return elem_[i][j]; }
   
    inline T* getElements() const { return elem_; } 
 
    void print() const {
      for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++)
          std::cout << elem_[i*N+j] << " ";
      }
    }

    /**
      * Sets the matrix element at position (i,j) to x. Indexes are checked.
      */
    void setElement(int i, int j, T x) {
      if (((i>=0) && (i < N)) && ((j>=0) && (j<N)))
        elem_[i*N+j] = x;
      else {
        std::cerr << "Invalid index in TinyMatrix: \n"
                  << "Matrix dimensions are ( " << N<< " x " << N << ")\n"
                  << "given indexes are [" << i << "][" << j << "].\n";
        boost::throw_exception(std::out_of_range("invalid index"));
      }
    }
  
    /**
      * Sets the matrix element at position (i,j) to x. Indexes are NOT checked.
      */
    inline void setElement_nc(int i, int j, T x) { elem_[i*N+j]=x; }
    
    /**
     * Computes v' * <this Matrix>
     */
    TinyVector<T,N> lmult(const TinyVector<T,N>& v) const {
      TinyVector<T,N> res;
      int index;
      for (int i=0; i<N; i++) {
        res[i] = 0;
        index = i;
        for (int j=0; j<N; j++) {          
          res[i] += v[j]*elem_[index];
          index += N;
        }
      }
      return res;
    }
       
    /**
     * Computes (v' * m * v) 
     */
    inline T vec_mat_vec(const TinyVector<T,N>& v) const;
    
    /**
     * Computes (v1' * m * v2) 
     */
    inline T vec_mat_vec(const TinyVector<T,N>& v1, const TinyVector<T,N>& v2)
               const; 
    
    /**
     * unary minus sign 
     */
    inline SquareMatrix<T,N> operator -() const;
    inline T operator[](unsigned int pos) const { return elem_[pos]; }

    /**
     * Division by a scalar
     */   
    SquareMatrix<T,N> operator / (T div)
    {
      if (div == 0) {
        std::cerr << "Cannot divide matrix elements by 0!\n";
        boost::throw_exception(std::runtime_error("division by zero"));
      }
      SquareMatrix<T,N> res;
      T inv = 1./div;
      for (int i=0; i<N*N; i++)
          res.elem_[i] = inv*elem_[i];
      return res;
    }

    SquareMatrix<T,N> operator * (T factor) {
      SquareMatrix<T,N> res;
      for (int i=0; i<N*N; i++)
        res.elem_[i] = factor*elem_[i];
      return res;
    }

    /**
     * Returns the determinant of the matrix.
     */
    T det() const;
    
    bool is_symmetric() {
      for (int i=0; i<N; i++)
        for (int j=i+1; j<N; j++)
          if (elem_[i*N+j] != elem_[j*N+i]) return false;
      return true;
    }

    friend std::ostream& operator << (std::ostream& os,
              SquareMatrix<T,N>& m){
      for (int i=0; i<N; i++) 
        os << m.elem_[i] << " ";
      return os;
    }
    
    friend std::ostream& operator >> (std::ostream& os, 
              SquareMatrix<T,N>& m) {
      for (int i=0; i<N; i++) 
        os  >> m.elem_[i];
      return os;
    }
    
  private:
    T elem_[N*N];
};

template<class T, int N>
void SquareMatrix<T,N>::setMatrix(std::string& str, alps::Parameters parms) {
  Helper<T> myHelper(str,parms);
  int elemCount = myHelper.elemcount();
  switch (elemCount) {
    case 1:
      // build a diagonal matrix with all the same elements on the diagonal
      {
      T theElement = myHelper.getNextElement();
      for (int i=0; i<N*N; i++)
        elem_[i] = 0.;
      for (int j=0; j<N*N; j+=N)
        elem_[j] = theElement;
      }
      break;
      
    case N:
      // build a diagnonal matrix, elements are given in the string
      for (int i=0; i<N*N; i++)
        elem_[i]=0.;
      for (int j=0; j<N*N; j+=N)
        elem_[j] = myHelper.getNextElement();
      break;
      
    case N*(N+1)/2:
      // build a symmetric matrix
      for (int i=0; i<N; i++)
        for (int j=0; j<=i; j++)
          elem_[i*N+j] = myHelper.getNextElement();
      // set the remaining elements
      for (int i=0; i<N; i++)
        for (int j=i+1; j<N; j++)
          elem_[i*N+j] = elem_[j*N+i];
      break;
      
    case N*N:
      // build a "complete" matrix, reading every element from the string
      for (int i=0; i<N*N; i++)
        elem_[i] = myHelper.getNextElement();
      break;
     
    default:
      std::cerr << "Wrong number of matrix elements:\n"
                << "found " << elemCount << " elements in the input string\n"
                << "only 1, " << N << ", " << N*(N+1)/2 << " or " << N*N
                <<" are valid numbers of elements for the matrix.\n";
      boost::throw_exception(std::runtime_error("invalid input"));
  }
};

template<class T, int N>
void SquareMatrix<T,N>::setMatrix(T x) {
  for (int i=0; i<N*N; i++)
    elem_[i] = x;
}

template <class T, int N>
struct squMat {
  // could be optimized ... 
  static T det(const SquareMatrix<T,N>& m) {
    SquareMatrix<T,N-1> minor;
    T det;
    for (int i=0; i<N-1; i++)
      for (int j=1; j<N; j++)
        minor.setElement_nc(i,j-1,m.getElement_nc(i,j));
    det = m.getElement_nc(N-1,0)*squMat<T,N-1>::det(minor);
    for (int s=0; s<N-1; s++) {
      for (int z=0; z<N-1; z++)
        minor.setElement_nc(z,s,m.getElement_nc(z,s));
      if (s & 1)
        det += m.getElement_nc(N-1,s+1)*squMat<T,N-1>::det(minor);
      else
        det -= m.getElement_nc(N-1,s+1)*squMat<T,N-1>::det(minor);
    }
    if (N & 1)
      return det;
    return -det;
  }

  static T vmv(const SquareMatrix<T,N>& m, const TinyVector<T,N>& v) {
    TinyVector<T,N> tmp = m.lmult(v);
    return dot(tmp,v);
  }
};

template <class T>
struct squMat<T,3> {
  static T det(const SquareMatrix<T,3>& m)
  { return (m[0]*m[4]*m[8] - m[0]*m[5]*m[7] + m[1]*m[5]*m[6] - m[1]*m[3]*m[8]
           +m[2]*m[3]*m[7] - m[2]*m[4]*m[6]); }
  static T vmv(const SquareMatrix<T,3>& m, const TinyVector<T,3>& v1,
               const TinyVector<T,3>& v2)
  { return (  v1[0] * ( m[0]*v2[0] + m[1]*v2[1] + m[2]*v2[2] )
            + v1[1] * ( m[3]*v2[0] + m[4]*v2[1] + m[5]*v2[2] )
            + v1[2] * ( m[6]*v2[0] + m[7]*v2[1] + m[8]*v2[2] ) ); }
};

template <class T>
struct squMat<T,2> {
  static T det(const SquareMatrix<T,2>& m) { return (m[0]*m[3] - m[1]*m[2]); }
  static T vmv(const SquareMatrix<T,2>& m, const TinyVector<T,2>& v1, 
               const TinyVector<T,2>& v2)
  { return (v1[0]*(m[0]*v2[0] + m[1]*v2[1]) + v1[1]*(m[2]*v2[0] + m[3]*v2[1]));}
};

template <class T>
struct squMat<T,1> {
  static T det(const SquareMatrix<T,1>& m) { return m[0]; }
  static T vmv(const SquareMatrix<T,1>& m, const TinyVector<T,1>& v1, 
               const SquareMatrix<T,1>& v2) { return (v1[0]*m[0]*v2[0]); }
};

template <class T>
struct squMat<T,0> {
  static T det(const SquareMatrix<T,0>&) { return 0; }
  static T vmv(const SquareMatrix<T,0>&, const TinyVector<T,0>&,
               const TinyVector<T,0>&) { return 0; }
};

template<class T, int N>
T det(const SquareMatrix<T,N>& m) { return squMat<T,N>::det(m); }

template<class T, int N>
T SquareMatrix<T,N>::det() const { return squMat<T,N>::det(*this); }

template<class T, int N>
T SquareMatrix<T,N>::vec_mat_vec(const TinyVector<T,N>& v) const
  { return squMat<T,N>::vmv(*this,v,v); }

template<class T, int N>
T SquareMatrix<T,N>::vec_mat_vec(const TinyVector<T,N>& v1,
               const TinyVector<T,N>& v2) const
  { return squMat<T,N>::vmv(*this,v1,v2); }

template<class T, int N>
inline SquareMatrix<T,N> SquareMatrix<T,N>::operator - () const {
  SquareMatrix<T,N> res;
  for(int i=0; i<N*N; i++) res.elem_[i] = - elem_[i];
  return res;
}
  
/**
 * Subclass for diagnonal Matrices
 */
template<class T, int N>
class DiagMatrix {
  public:
    /**
      * Constructer which sets all matrix elements to x.
      */
    DiagMatrix(T x=0.) { DiagMatrix<T,N>::setMatrix(x); }

    /** 
      * Constructer which parces the given string, evaluates the substrings and 
      * sets the matrix, in row-first manner. The substrings are delimited by 
      * white spaces.
      *
      * \param str the substring to be parced.
      * \param parms contains the values for expressions contained in the 
      *              string.
      */
    DiagMatrix(std::string& str,alps::Parameters parms) 
      { DiagMatrix<T,N>::setMatrix(str,parms); }

    static const int dim = N;
    static const bool allow_cluster_update = false;
 
    void setMatrix(std::string&, alps::Parameters); 
    void setMatrix(T x);

    /** 
     * Returns the element at position (i,j) of the matrix. The index bounds
     * are checked, an exception is thrown if the given indexes are not valid. 
     */
    T getElement(int i, int j) const {
      if (((i>=0) && (i < N)) && ((j>=0) && (j<N)))
        if (i==j)
          return elem_[i];
        else
          return 0.;
      else {
        std::cerr << "Invalid index in TinyMatrix: \n"
                  << "Matrix dimensions are ( " << N << " x " << N << ")\n"
                  << "given indexes are [" << i << "][" << j << "].\n";
        boost::throw_exception(std::out_of_range("invalid index"));
      }
    }
   
    /**
     * Returns the element at position (i,j) of the matrix. The index bounds are
     * NOT checked -- use this function only if you are sure that the indexes 
     * are ok!
     * The version is faster than getElement(int i, int j).
     */

    inline T getElement_nc(int i, int j) const { 
      if (i==j) return elem_[i];
      else return 0;
    }
   
    void print() const {
      for (int i=0; i<N; i++) 
        std::cout << elem_[i] << " ";
    }

    /**
      * Sets the matrix element at position (i,j) to x. Indexes are checked.
      */
    void setElement(int i, int j, T x) {
      if (((i>=0) && (i < N)) && ((j>=0) && (j<N)))
        if (i==j)
          elem_[i] = x;
        else {
          std::cerr << "Invalid index argument in DiagMatrix:\n"
                    << "cannot set off-diagonal element with indexes\n"
                    << "(" << i << "," << j << ") to value other than 0.\n";
          boost::throw_exception(std::runtime_error("invalid index arguments"));
        } 
      else {
        std::cerr << "Invalid index in TinyMatrix: \n"
                  << "Matrix dimensions are ( " << N << " x " << N << ")\n"
                  << "given indexes are [" << i << "][" << j << "].\n";
        boost::throw_exception(std::out_of_range("invalid index"));
      }
    }

    /**
      * Sets the matrix element at position (i,j) to x. Indexes are NOT checked.
      */
    inline void setElement_nc(int i, int j, T x) { 
      if (i==j) elem_[i]=x; 
    }

    /**
     * Computes (v' * m * v) 
     */
    inline T vec_mat_vec(const TinyVector<T,N>& v) const;
    
    /**
     * Computes (v1' * m * v2) 
     */
    inline T vec_mat_vec(const TinyVector<T,N>& v1, const TinyVector<T,N>& v2) 
               const;
     
    T det() const;
    
    bool is_symmetric() { return true; }

    /**
     * unary minus sign 
     */
    inline DiagMatrix<T,N> operator -() const;

    DiagMatrix<T,N> operator / (T div)
    {
      if (div == 0) {
        std::cerr << "Cannot divide matrix elements by 0!\n";
        boost::throw_exception(std::runtime_error("division by zero"));
      }
      DiagMatrix<T,N> res;
      T inv = 1./div;
      for (int i=0; i<N; i++)
        res.elem_[i] = inv*elem_[i];
      return res;
    }

    DiagMatrix<T,N> operator * (T factor) {
      DiagMatrix<T,N> res;
      for (int i=0; i<N; i++)
        res.elem_[i] = factor*elem_[i];
      return res;
    }


    friend std::ostream& operator << (std::ostream& os,
              const DiagMatrix<T,N> m){
      for (int i=0; i<N; i++) 
        os << m.elem_[i] << " ";
      return os;
    }

    friend std::ostream& operator >> (std::ostream& os, 
              DiagMatrix<T,N>& m) {
      for (int i=0; i<N; i++) {
        os >> m.elem_[i];
      }
      return os;
    }
    
    private:
      T elem_[N];
};

template<class T, int N>
void DiagMatrix<T,N>::setMatrix(std::string& str,alps::Parameters parms) {
  Helper<T> myHelper(str,parms);
  int elemCount = myHelper.elemcount();
  switch (elemCount) {
    case 1:
      // build a diagonal matrix with all the same elements on the diagonal
      {
      T theElement = myHelper.getNextElement();
      for (int i=0; i<N; i++)
        elem_[i] = theElement;
      }
      break;
      
    case N:
      // build a diagnonal matrix, elements are given in the string
      for (int i=0; i<N; i++)
        elem_[i] = myHelper.getNextElement();
      break;

    default:
      std::cerr << "Wrong number of matrix elements:\n"
                << "found " << elemCount << " elements in the input string\n"
                << "only 1, " << N << ", " << N*(N+1)/2 << " or " << N*N
                <<" are valid numbers of elements for the matrix.\n";
      boost::throw_exception(std::runtime_error("invalid input"));
  }
};

template<class T, int N>
void DiagMatrix<T,N>::setMatrix(T x) {
  for (int i=0; i<N; i++) {
    elem_[i] = x;
  }
}

template<class T, int N>
DiagMatrix<T,N> DiagMatrix<T,N>::operator - () const {
  DiagMatrix<T,N> res;
  for (int i=0; i<N; i++) 
    res.elem_[i] = -elem_[i];
  return res;
};    
 
template<class T, int N, int I>
struct diagMeta {
  static T vmv(const TinyVector<T,N>& v, const T elem[N]) {
    return (diagMeta<T,N,I-1>::vmv(v,elem) + v[I-1]*v[I-1]*elem[I-1]);
  }
  static T vmv(const TinyVector<T,N>& v1, const T elem[N],
               const TinyVector<T,N>& v2) {
    return (diagMeta<T,N,I-1>::vmv(v1,elem,v2)+v1[I-1]*elem[I-1]*v2[I-1]);
  }
  static T det(const T elem[N]) {
    return (diagMeta<T,N,I-1>::det(elem)*elem[I-1]);
  }
};

template <class T, int N>
struct diagMeta<T,N,0> {
  static T vmv(const TinyVector<T,N>& , const T[N]) { return 0.0; }
  static T vmv(const TinyVector<T,N>& , const T[N], 
               const TinyVector<T,N>& ) { return 0.0; }
  static T det(const T[N]) { return 1.; }    
};

template<class T, int N>
inline T DiagMatrix<T,N>::vec_mat_vec(const TinyVector<T,N>& v) const {
  return diagMeta<T,N,N>::vmv(v,elem_);
}    

template<class T, int N>
inline T DiagMatrix<T,N>::vec_mat_vec(const TinyVector<T,N>& v1,
                    const TinyVector<T,N>& v2) const {
  return diagMeta<T,N,N>::vmv(v1,elem_,v2);
}

template<class T, int N>
T DiagMatrix<T,N>::det() const { 
  return diagMeta<T,N,N>::det(elem_); 
}
    
/**
 * Subclass for Matrices with only one element 
 * (i.e. diagonal Matrices with all the same elements on the diagonal)
 */
template<class T, int N>
class MIdMatrix {
  public:
    /**
      * Constructer which sets all matrix elements to x.
      */
    MIdMatrix(T x=0.) { elem_ = x; }

    /** 
      * Constructer which parces the given string, evaluates the substrings and 
      * sets the matrix, in row-first manner. The substrings are delimited by 
      * white spaces.
      *
      * \param str the substring to be parced.
      * \param parms contains the values for expressions contained in the 
      *              string.
      */
    MIdMatrix(std::string& str,alps::Parameters parms) { setMatrix(str,parms); }
 
    static const int dim = N;  
    static const bool allow_cluster_update = true;
 
    inline void setMatrix(std::string&,alps::Parameters); 
    inline void setMatrix(T x) { elem_ = x; }

    /** 
     * Returns the element at position (i,j) of the matrix. The index bounds
     * are checked, an exception is thrown if the given indexes are not valid. 
     */
    T getElement(int i, int j) const {
      if (((i>=0) && (i < N)) && ((j>=0) && (j<N))) {
        if (i==j)
          return elem_;
        else
          return 0.;
      } else {
        std::cerr << "Invalid index in TinyMatrix: \n"
                  << "Matrix dimensions are ( " << N << " x " << N << ")\n"
                  << "given indexes are [" << i << "][" << j << "].\n";
        boost::throw_exception(std::out_of_range("invalid index"));
      }
    }

    /**
     * Returns the element at position (i,j) of the matrix. The index bounds are
     * NOT checked -- use this function only if you are sure that the indexes 
     * are ok!
     * The version is faster than getElement(int i, int j).
     */
    inline T getElement_nc(int , int ) const { return elem_; }
    
    void print() const {
      std::cout << elem_;
    }

    /**
      * Sets the matrix element at position (i,j) to x. Indexes are checked.
      */
    void setElement(int i, int j, T x) {
      if (((i>=0) && (i < N)) && ((j>=0) && (j<N)))
        if (i==j)
          elem_ = x;
        else {
          std::cerr << "Invalid index argument in DiagMatrix:\n"
                    << "cannot set off-diagonal element with indexes\n"
                    << "(" << i << "," << j << ") to value other than 0.\n";
          boost::throw_exception(std::runtime_error("invalid index arguments"));
        } 
      else {
        std::cerr << "Invalid index in TinyMatrix: \n"
                  << "Matrix dimensions are ( " << N << " x " << N << ")\n"
                  << "given indexes are [" << i << "][" << j << "].\n";
        boost::throw_exception(std::out_of_range("invalid index"));
      }
    }

  
    /**
      * Sets the matrix element at position (i,j) to x. Indexes are NOT checked.
      */
    inline void setElement_nc(int i, int j, T x) { elem_ = x; }

    /**
     * Computes (v' * m * v)
     */
    inline T vec_mat_vec(const TinyVector<T,N>& v) const 
    { return (elem_ * dot(v,v)); };
 
    /**
     * Computes (v1' * m * v2)
     */
    inline T vec_mat_vec(const TinyVector<T,N>& v1, const TinyVector<T,N>& v2) const
    { return elem_*dot(v1,v2); }
     
    T det() const;

    bool is_symmetric() { return true; }

    /**
     * unary minus sign 
     */
    inline MIdMatrix<T,N> operator -() const {
      MIdMatrix<T,N> res;
      res.elem_ = -elem_;
      return res;
    }

    MIdMatrix<T,N> operator / (T div)
    {
      if (div == 0) {
        std::cerr << "Cannot divide matrix elements by 0!\n";
        boost::throw_exception(std::runtime_error("division by zero"));
      }
      MIdMatrix<T,N> res;
      res.elem_ = elem_/div;
      return res;
    }
   
    MIdMatrix<T,N> operator * (T factor)
    {
      MIdMatrix<T,N> res;
      res.elem_ = elem_ * factor;
      return res;
    }

    friend std::ostream& operator << (std::ostream& os,
              const MIdMatrix<T,N> m){
      os << m.elem_ << " ";
      return os;
    }

    friend std::ostream& operator >> (std::ostream& os, 
              MIdMatrix<T,N>& m) {
      os >> m.elem_;
      return os;
    }
    
    private:
      T elem_;
};


template<class T, int N>
void MIdMatrix<T,N>::setMatrix(std::string& str,alps::Parameters parms) {
  Helper<T> myHelper(str,parms);
  int elemCount = myHelper.elemcount();
  if (elemCount == 1)
    elem_ = myHelper.getNextElement();
  else {
    std::cerr << "Wrong number of matrix elements:\n"
              << "found " << elemCount << " elements in the input string\n"
              << "only 1, " << N << ", " << N*(N+1)/2 << " or " << N*N
              << " are valid numbers of elements for the matrix.\n";
    boost::throw_exception(std::runtime_error("invalid input"));
  }
};

// modifications in the implementation of the const matrix, 
// which will hopefully become faster like this.
template<class T, int N, int I>
struct constMeta {
  static T len2(const TinyVector<T,N>& v)
    { return v[I-1]*v[I-1]+constMeta<T,N,I-1>::len2(v); }
  static T power(const T base)
    { return base*(constMeta<T,N,I-1>::power(base)); }
};

template<class T, int N>
struct constMeta<T,N,0> {
  static T len2(const TinyVector<T,N>&) { return 0.; }
  static T power(const T) { return 1.; }
};

template<class T, int N>
inline T MIdMatrix<T,N>::det() const 
  { return constMeta<T,N,N>::power(elem_); }  

#endif //ALPS_APPLICATIONS_MC_SPIN_MATRICES_H_
