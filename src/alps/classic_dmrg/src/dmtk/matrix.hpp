// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/matrix.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#ifndef __DMTK_MATRIX_H__
#define __DMTK_MATRIX_H__

// Generic Matrix class using templates
// Matrix is a vector subclass
// The data is stored in Fortran style: column0.column1....column(n-1) 
// for compatibility with numeric software like BLAS and LAPACK
// To retrieve the data in Fortran style use m(col,row)
// To retrieve the data in C style use m[row][col]

#include "dmtk/output.hpp"
#include <vector>
#include <iosfwd>
#include <iostream>
#include <utility>
#include "dmtk/conj.hpp"
#include "dmtk/slice_iter.hpp"
#include "dmtk/gslice_iter.hpp"
#include "dmtk/vector.hpp"
#include "dmtk/range.hpp"
#include "dmtk/array_util.hpp"

namespace dmtk
{
  using uint = unsigned int;

template<class T>
class Matrix : public std::vector<T>
{
  private:
    size_t num_rows = 0, num_cols = 0;

  public:
    using vector = std::vector<T>;
    using value_type = T;

    Matrix() = default;

    Matrix(size_t ncols, size_t nrows)
      : vector(nrows*ncols), num_rows(nrows), num_cols(ncols) {}
    Matrix(size_t ncols, size_t nrows, const T& v)
      : vector(nrows*ncols, v), num_rows(nrows), num_cols(ncols) {}
    Matrix(size_t ncols, size_t nrows, const T* v)
      : vector(nrows*ncols, v), num_rows(nrows), num_cols(ncols) {}

    Matrix(const Matrix<T>& m) = default;

    Matrix(cgslice_iter<T> s)
      : vector(s.size1()*s.size2()), num_rows(s.size2()), num_cols(s.size1())
      {
        for(size_t i = 0; i < rows(); i++)
          for(size_t j = 0; j < cols(); j++)
            this->operator()(j,i) = s(j,i);
      }

    Matrix(gslice_iter<T> s)
      : vector(s.size1()*s.size2()), num_rows(s.size2()), num_cols(s.size1())
      {
        for(size_t i = 0; i < rows(); i++)
          for(size_t j = 0; j < cols(); j++)
            this->operator()(j,i) = s(j,i);
      }

//  Operators

    T operator()(size_t col, size_t row) const; // Fortran style
    T& operator()(size_t col, size_t row);
    slice_iter<T> operator[] (size_t row); // C style
    cslice_iter<T> operator[] (size_t row) const; // C style

//  Asignment 
 
    Matrix& operator=(const Matrix<T> &);

    Matrix& operator=(cgslice_iter<T> s);
    Matrix& operator=(gslice_iter<T> s);
    Matrix& operator=(const T *v);
    Matrix& operator=(const T& v);

    Matrix& operator+=(const Matrix<T> &);
    Matrix& operator-=(const Matrix<T> &);
    Matrix& operator*=(const Matrix<T> &);
    Matrix& operator/=(const Matrix<T> &);

    Matrix& operator+=(const T &);
    Matrix& operator-=(const T &);
    Matrix& operator*=(const T &);
    Matrix& operator/=(const T &);

//  Iterators

    slice_iter<T> row(size_t);  // row iterator
    cslice_iter<T> row(size_t) const ;  // constant row iterator
    slice_iter<T> column(size_t);  // column iterator
    cslice_iter<T> column(size_t) const;  // constant column iterator
    T* array()             { return vector::data(); }
    const T* array() const { return vector::data(); }

//  Ranges

    slice_iter<T> operator()(Range, size_t);
    cslice_iter<T> operator()(Range, size_t) const;
    slice_iter<T> operator()(size_t, Range);
    cslice_iter<T> operator()(size_t, Range) const;
    gslice_iter<T> operator()(Range, Range);
    cgslice_iter<T> operator()(Range, Range) const;
    slice_iter<T> diagonal();
    cslice_iter<T> diagonal() const;

//  Methods
   
    Matrix<T>& randomize(int seed = -1)
    {
      // Fixed seed for ASLR-independent reproducibility (see the matching
      // comment on Vector<T>::randomize()). The original pointer->long cast
      // seed was UB and perturbed DMRG outputs at the 8th digit between
      // process runs.
      static long idum = seed == -1 ? 42 : seed;
      for(int i = 0; i < cols(); i++)
        for(int j = 0; j < rows(); j++)
          operator()(i,j) = quickran(idum);
      return *this;
    }

    inline T trace() const 
    { 
      T r = T(0);
      for(uint i = 0; i < std::min(num_rows, num_cols); i++) 
        r += operator()(i,i);
      return r;
    }

    inline size_t rows() const { return num_rows; }
    inline size_t cols() const { return num_cols; }

//  Memory management

    void swap(Matrix& other) noexcept
    {
      vector::swap(other);
      std::swap(num_rows, other.num_rows);
      std::swap(num_cols, other.num_cols);
    }

    inline size_t memory() const { return this->capacity() * sizeof(T); }
    inline size_t usage() const  { return num_rows * num_cols * sizeof(T); }
    Matrix<T>& reshape(size_t ncols, size_t nrows);
    Matrix<T>& resize(size_t ncols, size_t nrows);
    Matrix<T>& resize(size_t ncols, size_t nrows, const T &);
    Matrix<T>& reserve(size_t ncols, size_t nrows);

    inline size_t size1() const { return num_cols; }
    inline size_t size2() const { return num_rows; }

//  Streams

    void read(std::istream& s);
    void write(std::ostream& s) const;
};

#include "dmtk/matrix_implement.hpp"

//////////////////////////////////////////////////////////

// Assignment

template<class T>
inline Matrix<T>& 
Matrix<T>::operator=(const Matrix<T> &m) 
{
  if(rows() != m.num_rows || cols() != m.num_cols)
    reshape(m.num_cols, m.num_rows);

  const T* orig = m.array();
  T *dest = array();

  size_t n = num_rows*num_cols;
  array_copy(n, orig, dest);
//  while(n--) *dest++ = *orig++; 

  return *this;
}

template<class T>
inline Matrix<T>& 
Matrix<T>::operator=(const T *v) 
{
  T *dest = array();
  size_t n = num_rows*num_cols;
  array_copy(n, v, dest);
//  while(n--) *dest++ = *v++; 

  return *this;
}

template<class T>
inline Matrix<T>& 
Matrix<T>::operator=(const T& v) 
{
  T *dest = array();
  size_t n = num_rows*num_cols;
  while(n--) *dest++ = v; 
  return *this;
}

template<class T>
inline Matrix<T>&
Matrix<T>::operator=(cgslice_iter<T> mm)
{
  reshape(mm.size1(), mm.size2());
  for(size_t i = 0; i < rows(); i++)
    for(size_t j = 0; j < cols(); j++)
      this->operator()(j,i) = mm(j,i);

  return *this;
}

template<class T>
inline Matrix<T>&
Matrix<T>::operator=(gslice_iter<T> mm)
{
  reshape(mm.size1(), mm.size2());
  for(size_t i = 0; i < rows(); i++)
    for(size_t j = 0; j < cols(); j++)
      this->operator()(j,i) = mm(j,i);

  return *this;
}

//////////////////////////////////////////////////////////

#define BINARY_OP(op,ap) \
template<class T> \
inline Matrix<T>& \
Matrix<T>::op(const Matrix<T>& mm) \
{ \
  reshape(mm.cols(),mm.rows()); \
  const T* orig = mm.array(); \
  T* dest = array(); \
  int n = mm.cols()*mm.rows(); \
  while(n--) *dest++ ap *orig++; \
  return *this; \
} 

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

#define BINARY_OP(op,ap) \
template<class T> \
inline Matrix<T>& \
Matrix<T>::op(const T& v) \
{ \
  for(size_t i = 0; i < rows(); i++) \
    for(size_t j = 0; j < cols(); j++) \
      operator()(j,i) ap v; \
 \
  return *this; \
} 

BINARY_OP(operator+=,+=);
BINARY_OP(operator-=,-=);
BINARY_OP(operator/=,/=);
BINARY_OP(operator*=,*=);
#undef BINARY_OP

//////////////////////////////////////////////////////////

// Reference

template<class T>
inline slice_iter<T> 
Matrix<T>::row(size_t i) 
{
  return slice_iter<T>(this, std::slice(i, num_cols, num_rows));
}

template<class T>
inline cslice_iter<T> 
Matrix<T>::row(size_t i) const 
{
  return cslice_iter<T>(this, std::slice(i, num_cols, num_rows));
}

template<class T>
inline slice_iter<T> 
Matrix<T>::column(size_t i) 
{
  return slice_iter<T>(this, std::slice(i*num_rows, num_rows, 1));
}

template<class T>
inline cslice_iter<T> 
Matrix<T>::column(size_t i) const
{
  return cslice_iter<T>(this, std::slice(i*num_rows, num_rows, 1));
}

template<class T>
inline T 
Matrix<T>::operator()(size_t ncol, size_t nrow) const 
{
  return vector::operator[](ncol*num_rows+nrow);
}

template<class T>
inline T&
Matrix<T>::operator()(size_t ncol, size_t nrow) 
{
  return vector::operator[](ncol*num_rows+nrow);
}

template<class T>
inline slice_iter<T> 
Matrix<T>::operator[](size_t nrow) 
{
  return row(nrow);
}

template<class T>
inline cslice_iter<T> 
Matrix<T>::operator[](size_t nrow) const 
{
  return row(nrow);
}

//////////////////////////////////////////////////////////

// Memory management

template<class T>
Matrix<T>&
Matrix<T>::resize(size_t ncols, size_t nrows)
{
  if(nrows == num_rows && ncols == num_cols) return *this;

  Matrix<T> aux(*this);
  const size_t new_size = nrows * ncols;

  if (new_size > this->size())
    vector::resize(new_size);

  T *orig, *dest;
  for(uint col = 0 ; col < std::min(num_cols, ncols); col++){
    orig = &aux.column(col)[0];
    dest = &column(col)[0];
    for(uint row = 0 ; row < std::min(num_rows, nrows) ; row++){
        *dest++ = *orig++;
    }
  }
  num_cols = ncols;
  num_rows = nrows;

  return *this;
}

template<class T>
Matrix<T>&
Matrix<T>::resize(size_t ncols, size_t nrows, const T& v)
{
  if(nrows == num_rows && ncols == num_cols) return *this;

  Matrix<T> aux(*this);
  const size_t new_size = nrows * ncols;

  if (new_size > this->size())
    vector::resize(new_size, v);

  T *orig, *dest;
  for(uint col = 0 ; col < std::min(num_cols, ncols); col++){
    orig = &aux.column(col)[0];
    dest = &column(col)[0];
    for(uint row = 0 ; row < std::min(num_rows, nrows) ; row++){
        *dest++ = *orig++;
    }
  }
  num_cols = ncols;
  num_rows = nrows;

  return *this;
}

template<class T>
Matrix<T>&
Matrix<T>::reshape(size_t ncols, size_t nrows)
{
  if(nrows == num_rows && ncols == num_cols) return *this;

  const size_t new_size = nrows * ncols;
  if (new_size > this->size())
    vector::resize(new_size);

  num_cols = ncols;
  num_rows = nrows;

  return *this;
}

template<class T>
Matrix<T>&
Matrix<T>::reserve(size_t ncols, size_t nrows)
{
  const size_t new_size = nrows * ncols;
  if (new_size > this->size()) {
    Matrix<T> aux(*this);
    vector::resize(new_size);

    T *orig, *dest;
    for(uint col = 0 ; col < num_cols; col++){
      orig = &aux.column(col)[0];
      dest = &column(col)[0];
      for(uint row = 0 ; row < num_rows ; row++){
          *dest++ = *orig++;
      }
    }
  }

  return *this;
}


////////////////////////////////////////////////////////
template<class T>
inline T
dot_product(const Matrix<T> &a, const Matrix<T> &b)
{
  using namespace std;
                                                                                
  if(a.size() != b.size())
     dmtk::err() << "** Warning: Matrix sizes do not comform\n";
                                                                                
  const T* pa = a.array();
  const T* pb = b.array();
                               
  int n = std::min(a.size(), b.size());
  return dot_product(n, pa, 1, pb, 1);
}
///////////////////////////////////////////////////////////////////////////

template<class T>
void
Matrix<T>::write(std::ostream &s) const
{
  s.write((const char *)&num_rows, sizeof(size_t));
  s.write((const char *)&num_cols, sizeof(size_t));

  for(uint col = 0; col < cols(); col++) 
    for(uint row = 0; row < rows(); row++){ 
      T v = operator()(col,row);
      s.write((const char *)&v, sizeof(T));
    }
}

template<class T>
void
Matrix<T>::read(std::istream &s)
{
  size_t nr, nc;
  s.read((char *)&nr, sizeof(size_t));
  s.read((char *)&nc, sizeof(size_t));
  resize(nc,nr);

  for(uint col = 0; col < cols(); col++) 
    for(uint row = 0; row < rows(); row++){ 
      T& v = operator()(col,row);
      s.read((char *)&v, sizeof(T));
    }
}

} // namespace dmtk 

#endif // __DMTK_MATRIX_H__ 
