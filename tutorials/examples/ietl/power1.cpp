/***************************************************************************
 * $Id: power1.cpp,v 1.5 2004/06/29 08:31:02 troyer Exp $
 *
 * A simple example of power method
 *
 * Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>
 *                            Matthias Troyer <troyer@comp-phys.org>
 *
* SPDX-License-Identifier: MIT
*
 *
 **************************************************************************/

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <ietl/interface/ublas.h>
#include <ietl/power.h>
#include <ietl/vectorspace.h>
#include <ietl/iteration.h>
#include <boost/random.hpp>
#include <boost/limits.hpp>


typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower> Matrix; 
typedef boost::numeric::ublas::vector<double> Vector;

int main () {
  int N = 10;
  Matrix mat(N, N);
  int n = 1;
  for(int i=0;i<N;i++)
    for(int j=0;j<=i;j++)
      mat(i,j) = n++;   
  
  std::cout << "Matrix: " << mat << std::endl;
  
  ietl::vectorspace<Vector> vec(N);
  boost::lagged_fibonacci607 gen;  
  int max_iter = std::numeric_limits<int>::max();
  double rel_tol = 5.*std::numeric_limits<double>::epsilon();
  double abs_tol = std::numeric_limits<double>::epsilon();
  ietl::basic_iteration<double> iter(max_iter, rel_tol, abs_tol);
  
  std::pair<double,Vector> result = ietl::power(mat, gen, iter, vec); 
  std::cout.precision(20);
  std::cout << "Eigenvalue: "<< result.first << std::endl;
  std::cout << "Eigenvector: " << result.second << std::endl;  
  return 0;
}
