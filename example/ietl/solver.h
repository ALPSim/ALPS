/***************************************************************************
 * $Id: solver.h,v 1.4 2004/03/04 12:34:58 troyer Exp $
 *
 * a LAPACK linear equation solver wrapper 
 *
 * Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>
 *                            Matthias Troyer <troyer@comp-phys.org>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
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

