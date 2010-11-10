/***************************************************************************
 * $Id: solver.h,v 1.4 2004/03/04 12:34:58 troyer Exp $
 *
 * a LAPACK linear equation solver wrapper 
 *
 * Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>
 *                            Matthias Troyer <troyer@comp-phys.org>
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

