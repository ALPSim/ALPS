/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_LAPACK_H
#define LOOPER_LAPACK_H

#include <cassert>
#include <cmath>
#include <complex>
#include <stdexcept>

#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/lapack/driver/syev.hpp>
#include <boost/numeric/bindings/lapack/driver/heev.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace looper {

//
// diagonalize matrix
//

template <class T, class R, class A, class Vector>
inline void diagonalize(
  boost::numeric::ublas::matrix<T,R,A>& a,
  Vector& w, bool need_eigenvectors = true)
{
  using namespace boost::numeric;

  BOOST_STATIC_ASSERT((boost::is_same<typename R::orientation_category,
    ublas::column_major_tag>::value));

  const char jobz = (need_eigenvectors ? 'V' : 'N');
  const char uplo = 'L';

  // call dispatcher
  int info = bindings::lapack::syev(jobz, uplo, a, w,
    bindings::lapack::optimal_workspace());
  if (info != 0) throw std::runtime_error("failed in heev");
}

template <class T, class R, class A, class Vector>
inline void diagonalize(
  boost::numeric::ublas::matrix<std::complex<T>,R,A>& a,
  Vector& w, bool need_eigenvectors = true)
{
  using namespace boost::numeric;

  BOOST_STATIC_ASSERT((boost::is_same<typename R::orientation_category,
    ublas::column_major_tag>::value));

  const char jobz = (need_eigenvectors ? 'V' : 'N');
  const char uplo = 'L';

  // call dispatcher
  int info = bindings::lapack::heev(jobz, uplo, a, w,
    bindings::lapack::optimal_workspace());
  if (info != 0) throw std::runtime_error("failed in heev");
}


//
// solve_llsp: solving linear least-squares problem
//

template <typename T, typename R, typename A>
inline double solve_llsp(
  const boost::numeric::ublas::matrix<T, R, A>& a,
  const boost::numeric::ublas::vector<T>& b,
  boost::numeric::ublas::vector<T>& x,
  double tol = 1.0e-10)
{
  using namespace boost::numeric;
#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
  using ublas::norm_inf; using ublas::norm_2; using ublas::prod;
#endif

  BOOST_STATIC_ASSERT((boost::is_same<typename R::orientation_category,
    ublas::column_major_tag>::value));

  typedef T value_type;
  typedef ublas::matrix<value_type, R, A> matrix_type;
  typedef ublas::vector<value_type> vector_type;

  int const m = bindings::size_row(a);
  int const n = bindings::size_column(a);
  int const min_mn = std::min(m, n);

  // temporary storage
  matrix_type at; at = a;
  matrix_type u(m, min_mn);
  matrix_type vt(min_mn, n);
  vector_type s(min_mn);

  // call SVD
  int info = bindings::lapack::gesvd('S','S',at, s, u, vt);
  if (info != 0) throw std::runtime_error("failed in gesvd");

  // inverse S
  double smax_inv = 1.0 / norm_inf(s);
  for (int i = 0; i < min_mn; ++i)
    s(i) = (smax_inv * s(i) > tol) ? (1.0 / s(i)) : 0.0;

  for (int i = 0; i < n; ++i) {
    x(i) = 0.0;
    for (int j = 0; j < min_mn; ++j)
      for (int k = 0; k < m; ++k) x(i) += vt(j,i) * s(j) * u(k,j) * b(k);
  }

  // return residual
  return norm_2(prod(a,x) - b);
}

} // end namespace looper

#endif // ! LOOPER_LAPACK_H
