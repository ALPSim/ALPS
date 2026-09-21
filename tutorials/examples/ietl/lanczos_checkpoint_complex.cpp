/***************************************************************************
 * $Id: lanczos1.cpp,v 1.9 2004/06/29 08:31:02 troyer Exp $
 *
 * An example of the Lanczos method for the calculation of n lowest eigenvalues.
 *
 * Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>
 *                            Matthias Troyer <troyer@comp-phys.org>
 *
* SPDX-License-Identifier: MIT
*
 *
 **************************************************************************/

#include <alps/osiris.h>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <ietl/interface/ublas.h>
#include <ietl/vectorspace.h>
#include <ietl/lanczos.h>
#include <boost/random.hpp>
#include <boost/limits.hpp>
#include <cmath>
#include <limits>
#include <complex>
#include <alps/osiris/boost/ublas.h>

typedef boost::numeric::ublas::hermitian_matrix<std::complex<double> > Matrix; 
typedef boost::numeric::ublas::vector<std::complex<double> > Vector;

int main() {
  // Creation of an example matrix:
  int N = 10;
  Matrix mat(N, N);
  int n = 1;
  for(int i=0;i<N;i++)
    for(int j=0;j<=i;j++)
      mat(i,j) = n++;    
  std::cout << "\n" << "Printing matrix\n";
  std::cout << "--------------------------------\n\n";
  std::cout << mat << std::endl;

  typedef ietl::vectorspace<Vector> Vecspace;
  typedef boost::lagged_fibonacci607 Gen;  

  Vecspace vec(N);
  Gen mygen;
  ietl::lanczos<Matrix,Vecspace> lanczos(mat,vec);

  // Creation of an iteration object:  
  int max_iter = 10*N;  
  double rel_tol = 500*std::numeric_limits<double>::epsilon();
  double abs_tol = std::pow(std::numeric_limits<double>::epsilon(),2./3);  
  std::cout << "Computation of 2 lowest converged eigenvalues\n\n";
  std::cout << "-----------------------------------\n\n";
  int n_lowest_eigenval = 2;
  std::vector<double> eigen;
  std::vector<double> err;
  std::vector<int> multiplicity;  
  ietl::lanczos_iteration_nlowest<double> 
    iter(max_iter, n_lowest_eigenval, rel_tol, abs_tol);
  try{
    lanczos.calculate_eigenvalues(iter,mygen);
    eigen = lanczos.eigenvalues();
    err = lanczos.errors();
    multiplicity = lanczos.multiplicities();
    std::cout<<"number of iterations: "<<iter.iterations()<<"\n";
  }
  catch (std::runtime_error& e) {
    std::cout << e.what() << "\n";
  } 
  
  // Printing eigenvalues with error & multiplicities:  
  std::cout << "#        eigenvalue            error         multiplicity\n";  
  std::cout.precision(10);
  for (unsigned int i=0;i<eigen.size();++i) 
    std::cout << i << "\t" << eigen[i] << "\t" << err[i] << "\t" 
          << multiplicity[i] << "\n";

  alps::OXDRFileDump ar ("lancos_complex.dump",0);
  ar << lanczos;
  return 0;
}
