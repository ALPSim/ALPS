// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/matrix_implement.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#ifndef __DMTK_MATRIX_IMPLEMENT_H__
#define __DMTK_MATRIX_IMPLEMENT_H__

// Matrix member impls that need slice_iter / gslice_iter types fully
// defined. Historically this file also held the expression-template
// MATRIX_* / SCALAR_MATRIX_OP macros, MatrixVectorProduct, MatrixTranspose,
// tensor(), I<>(), inverse(), etc. — all dead once meta.h was retired.

/////////////////////////////////////////////////////////////////

template<class T>
slice_iter<T>
Matrix<T>::operator()(Range col_range, size_t row)
{
  std::slice s(col_range.start()*num_rows + row,
          col_range.size(),
          col_range.stride()*num_rows);

  slice_iter<T> iter(this, s);
  return iter;
}

template<class T>
cslice_iter<T>
Matrix<T>::operator()(Range col_range, size_t row) const
{
  std::slice s(col_range.start()*num_rows + row,
          col_range.size(),
          col_range.stride()*num_rows);

  cslice_iter<T> iter(this, s);
  return iter;
}

template<class T>
slice_iter<T>
Matrix<T>::operator()(size_t col, Range row_range)
{
  std::slice s(col*num_rows+row_range.start(),
          row_range.size(),
          row_range.stride());

  slice_iter<T> iter(this, s);
  return iter;
}

template<class T>
cslice_iter<T>
Matrix<T>::operator()(size_t col, Range row_range) const
{
  std::slice s(col*num_rows+row_range.start(),
          row_range.size(),
          row_range.stride());

  cslice_iter<T> iter(this, s);
  return iter;
}

template<class T>
gslice_iter<T>
Matrix<T>::operator()(Range col_range, Range row_range)
{
  std::slice s_col(col_range.start()*num_rows,
           col_range.size(),
           col_range.stride()*num_rows);
  std::slice s_row(row_range.start(),
           row_range.size(),
           row_range.stride());
  gslice_iter<T> iter(this, s_col, s_row);
  return iter;
}

template<class T>
cgslice_iter<T>
Matrix<T>::operator()(Range col_range, Range row_range) const
{
  std::slice s_col(col_range.start()*num_rows,
           col_range.size(),
           col_range.stride()*num_rows);
  std::slice s_row(row_range.start(),
           row_range.size(),
           row_range.stride());
  cgslice_iter<T> iter(this, s_col, s_row);
  return iter;
}

template<class T>
slice_iter<T>
Matrix<T>::diagonal()
{
  std::slice s(0, std::min(num_rows,num_cols), num_rows+1);

  slice_iter<T> iter(this, s);
  return iter;
}

template<class T>
cslice_iter<T>
Matrix<T>::diagonal() const
{
  std::slice s(0, std::min(num_rows,num_cols), num_rows+1);

  cslice_iter<T> iter(this, s);
  return iter;
}

template<class T>
slice_iter<T>
inline diagonal(Matrix<T> &m)
{
  return m.diagonal();
}

template<class T>
inline T trace(Matrix<T> &m)
{
  return m.trace();
}

// Elementwise conjugate — returns a fresh Matrix<T>. Previously this
// was an expression-template helper defined in matrix_implement.hpp via
// MATRIX_METHOD(conj, MatrixConj); it is now materialised directly so
// that `dmtk::product(..., hc=true)` keeps compiling.
template<class T>
inline Matrix<T> conj(const Matrix<T>& m)
{
  Matrix<T> r(m.cols(), m.rows());
  for (size_t j = 0; j < m.cols(); ++j)
    for (size_t i = 0; i < m.rows(); ++i)
      r(j, i) = conj(m(j, i));
  return r;
}

#endif /* __DMTK_MATRIX_IMPLEMENT_H__ */
