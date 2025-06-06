/***************************************************************************
 * $Id: lanczos1.cpp,v 1.9 2004/06/29 08:31:02 troyer Exp $
 *
 * An example of the Lanczos method for the calculation of n lowest eigenvalues.
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
 *
 **************************************************************************/

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <ietl/interface/ublas.h>
#include <ietl/vectorspace.h>
#include <ietl/lanczos.h>
#include <boost/random.hpp>
#include <boost/limits.hpp>
#include <cmath>
#include <limits>

typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::lower> Matrix; 
typedef boost::numeric::ublas::vector<double> Vector;

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

  // Another set of computations for more converged eigenvalues: 
  std::cout << "\nMore converged eigenvalues\n\n";
  std::cout << "-----------------------------------\n\n";
  n_lowest_eigenval = 7;  
  try {
    ietl::lanczos_iteration_nlowest<double> 
      iter(max_iter, n_lowest_eigenval, rel_tol, abs_tol);
    lanczos.more_eigenvalues(iter); 
    eigen = lanczos.eigenvalues();
    err = lanczos.errors();
    multiplicity = lanczos.multiplicities();
    std::cout<<"number of iterations: "<<iter.iterations()<<"\n";
  }
  catch (std::runtime_error& e) {
    std::cout << e.what() << "\n";
  } 
  
  std::cout << "#        eigenvalue            error         multiplicity\n";
  for (unsigned int i=0;i<eigen.size();++i) 
    std::cout << i << "\t" << eigen[i] << "\t" << err[i] << "\t" 
          << multiplicity[i] << "\n";  

  // call of eigenvectors function follows:   
  std::cout << "\nEigen vectors computations for 3 lowest eigenvalues:\n\n";  
  std::vector<double>::iterator start = eigen.begin();  
  std::vector<double>::iterator end = eigen.begin()+3;
  std::vector<Vector> eigenvectors; // for storing the eigen vectors. 
  ietl::Info<double> info; // (m1, m2, ma, eigenvalue, residualm, status).
  
  try {
    lanczos.eigenvectors(start,end,std::back_inserter(eigenvectors),info,mygen); 
  }
  catch (std::runtime_error& e) {
    std::cout << e.what() << "\n";
  }  
  
  std::cout << "Printing eigenvectors:\n\n"; 
  for(std::vector<Vector>::iterator it = eigenvectors.begin();it!=eigenvectors.end();it++){
    std::copy((it)->begin(),(it)->end(),std::ostream_iterator<double>(std::cout,"\n"));
    std::cout << "\n\n";
  }
  std::cout << " Information about the eigenvector computations:\n\n";
  for(int i = 0; i < info.size(); i++) {
    std::cout << " m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
          << info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
          << i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
          << info.residual(i) << " error_info(" << i+1 << "): "
          << info.error_info(i) << "\n\n";
  }
  return 0;
}
