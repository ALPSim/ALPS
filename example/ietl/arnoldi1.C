/***************************************************************************
 * $Id: arnoldi1.cpp,v 1.9 2004/06/29 08:31:02 troyer Exp $
 *
 * An example of the Arnoldi iteration for non-Hermitian matrices
 *
 * Copyright (C) 2001-2003 by Bela Bauer <bauerb@phys.ethz.ch>
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

#include <iostream>
#include <complex>

#include <boost/random.hpp>
#include <boost/limits.hpp>

#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/lapack/driver/geev.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

#include <boost/lambda/bind.hpp>
#include <boost/lambda/construct.hpp>

#include <ietl/interface/ublas.h>
#include <ietl/simple_arnoldi.h>
#include <ietl/vectorspace.h>

typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> Matrix; 
typedef boost::numeric::ublas::vector<double> Vector;

bool cmp(std::complex<double> a, std::complex<double> b)
{
	return real(a) > real(b);
}

int main () {
	boost::lagged_fibonacci607 gen;

   int N = 100;
   Matrix mat(N, N);
   int n = 1;
   for(int i=0;i<N;i++)
       for(int j=0;j<N;j++)
			mat(i,j) = gen(); //n++;
   
   ietl::vectorspace<Vector> vs(N);
   
   unsigned int nvals = 5;
   ietl::arnoldi_iteration<double> iter(50, nvals, 1e-8, 1e-8);
   
   ietl::simple_arnoldi<Matrix, ietl::vectorspace<Vector>, boost::lagged_fibonacci607> arni(mat, vs, gen);
   arni.calculate_eigenvalues(iter, true);
   
   for (int i = 0; i < nvals; ++i)
       std::cout << "Eigenvalue " << i << ": " << arni.get_eigenvalue(i) << std::endl;
   
    // Check with dense LAPACK
    {
     Matrix null1(N,N), null2(N,N);
     std::vector<std::complex<double> > evals(N);
     std::vector<double> ev1(N), ev2(N);
     boost::numeric::bindings::lapack::geev('N', 'N', mat, ev1, ev2, null1, null2);
     std::transform(ev1.begin(), ev1.end(), ev2.begin(), evals.begin(),
         boost::lambda::bind(boost::lambda::constructor<std::complex<double> >(), boost::lambda::_1, boost::lambda::_2));
         
     std::sort(evals.begin(), evals.end(), cmp);
     std::copy(evals.begin(), evals.begin()+nvals, std::ostream_iterator<std::complex<double> >(std::cout, " "));
     std::cout << std::endl;
    }

   return 0;
}
