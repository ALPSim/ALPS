/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_LAPACK_H
#define LOOPER_LAPACK_H

#include <cmath>
#include <complex>
#include <stdexcept>

#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/lapack/driver/syev.hpp>
#include <boost/numeric/bindings/lapack/driver/heev.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas.hpp>

//
// solve_llsp: solving linear least-squares problem
//

template <typename T, typename R, typename A>
inline T solve_llsp(
  const boost::numeric::ublas::matrix<T, R, A>& a,
  const boost::numeric::ublas::vector<T>& b,
  boost::numeric::ublas::vector<T>& x,
  double tol = 1.0e-10)
{
  using namespace boost::numeric;
#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
  using ublas::norm_inf; using ublas::norm_2; using ublas::prod;
#endif

  typedef T value_type;
  typedef ublas::matrix<value_type, R, A> matrix_type;
  typedef ublas::vector<value_type> vector_type;

  int const m = bindings::size_row(a);
  int const n = bindings::size_column(a);
  int const min_mn = std::min BOOST_PREVENT_MACRO_SUBSTITUTION (m, n);

  // temporary storage
  matrix_type at(a);
  matrix_type u(m, min_mn);
  matrix_type vt(min_mn, n);
  vector_type s(min_mn);

  // call SVD
  int info = bindings::lapack::gesvd('S','S',at, s, u, vt);
  if (info != 0) throw std::runtime_error("failed in gesvd");

  // compute the condition number to return lateron
  T cond = s(0) / s(min_mn-1);

  // inverse S
  double smax_inv = 1.0 / norm_inf(s);
  for (int i = 0; i < min_mn; ++i)
    s(i) = (smax_inv * s(i) > tol) ? (1.0 / s(i)) : 0.0;

  for (int i = 0; i < n; ++i) {
    x(i) = 0.0;
    for (int j = 0; j < min_mn; ++j)
      for (int k = 0; k < m; ++k) x(i) += vt(j,i) * s(j) * u(k,j) * b(k);
  }

  // return the condition number of the Matrix A, i.e. sigma(max)/sigma(min)
  return cond;
}

#endif // ! LOOPER_LAPACK_H
