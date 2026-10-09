/***************************************************************************
 * $Id: solver.h,v 1.4 2004/03/04 12:34:58 troyer Exp $
 *
 * a LAPACK linear equation solver wrapper 
 *
 * Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>
 *                            Matthias Troyer <troyer@comp-phys.org>
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
 **************************************************************************/

#include <ietl/traits.h>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/lapack/driver/sysv.hpp>
#include <boost/numeric/bindings/lapack/driver/hesv.hpp>
#include <complex>


template <class T> struct solver_helper 
{
  template <class M, class V>
  static void solve(M& m, V& v) { 
    boost::numeric::bindings::lapack::sysv('U',m,v);
  }
};

template <class T> struct solver_helper<std::complex<T> >
{
  template <class M, class V>
  static void solve(M& m, V& v) { 
    boost::numeric::bindings::lapack::hesv('U',m,v);
  }
};



template <class MATRIX, class VECTOR>
struct Solver
{
  typedef VECTOR vector_type;
  typedef typename vector_type::value_type scalar_type;
  typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
  typedef MATRIX matrix_type;
    
  void operator() (const matrix_type& mat, magnitude_type rho, const vector_type& x, vector_type& y) const {
    ietl::copy(x,y);
    matrix_type mat_ = mat -rho*boost::numeric::ublas::identity_matrix<scalar_type>(mat.size1());
    solver_helper<scalar_type>::solve(mat_,y);
  }
};

