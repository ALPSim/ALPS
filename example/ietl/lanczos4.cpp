/***************************************************************************
 * $Id: lanczos4.cpp,v 1.3 2004/06/29 09:27:48 troyer Exp $
 *
 * An example of the Lanczos method with fixed size T matrix.
 *
 * Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>
 *                            Matthias Troyer <troyer@comp-phys.org>
 *
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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
 *
 **************************************************************************/

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <ietl/interface/ublas.h>
#include <ietl/vectorspace.h>
#include <ietl/iteration.h>
#include <ietl/lanczos.h>
#include <boost/random.hpp>
#include <boost/limits.hpp>
#include <cmath>
#include <limits>

typedef std::complex<double> value_type;
typedef boost::numeric::ublas::hermitian_matrix<value_type, boost::numeric::ublas::lower> Matrix; 
typedef boost::numeric::ublas::vector<value_type> Vector;

int main() { 
  // Creation of a sample matrix:
  int N = 10; 
  Matrix mat(N, N);
  int n = 1;
  for(int i=0;i<N;i++)
    for(int j=0;j<=i;j++)
      mat(i,j) = n++;    
  std::cout << std::endl << "Printing matrix\n";
  std::cout << "--------------------------------\n\n";
  std::cout << mat << std::endl;
  
  typedef ietl::vectorspace<Vector> Vecspace;
  typedef boost::lagged_fibonacci607 Gen;  
  Vecspace vec(N);
  Gen mygen;
  ietl::lanczos<Matrix,Vecspace> lanczos(mat,vec);
  
  int max_iter = 2*N;  
  std::cout << "Computation of eigenvalues with fixed size of T-matrix\n\n";
  std::cout << "-----------------------------------\n\n";

  std::vector<double> eigen;
  std::vector<double> err;
  std::vector<int> multiplicity;
  ietl::fixed_lanczos_iteration<double> iter(max_iter);  
  try {
    lanczos.calculate_eigenvalues(iter,mygen);
    eigen = lanczos.eigenvalues();
    err = lanczos.errors();
    multiplicity = lanczos.multiplicities();
  }
  catch (std::runtime_error& e) {
    std::cout << e.what() << std::endl;
  } 
 // Printing eigenvalues with error & multiplicities:  
  std::cout << "#         eigenvalue         error         multiplicity\n";  
  std::cout.precision(10);
  for (unsigned int i=0;i<eigen.size();++i)
    std::cout << i << "\t" << eigen[i] << "\t" << err[i] << "\t" 
          << multiplicity[i] << "\n";
  
  // call of eigenvectors function follows:   
  std::cout << "\nEigen vectors computations for 3 lowest eigenvalues:\n\n";  
  std::vector<double>::iterator start = eigen.begin();
  std::vector<double>::iterator end = eigen.begin()+3;
  std::vector<Vector> eigenvectors; // for storing the eigen vector. 
  ietl::Info<double> info; // (m1, m2, ma, eigenvalue, residual, status).
  
  try {
    lanczos.eigenvectors(start, end, std::back_inserter(eigenvectors),info,mygen); 
  }
  catch (std::runtime_error& e) {
    std::cout << e.what() << std::endl;
  }  
  
  std::cout << "Printing eigen Vectors:" << std::endl << std::endl; 
  for(std::vector<Vector>::iterator it = eigenvectors.begin(); it != eigenvectors.end(); it++){
    std::copy((it)->begin(),(it)->end(),std::ostream_iterator<value_type>(std::cout,"\n"));
    std::cout << "\n\n";
  }
  std::cout << " Information about the eigen vectors computations:\n\n";
  for(int i = 0; i < info.size(); i++) {
    std::cout << " m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
          << info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
          << i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
          << info.residual(i) << " error_info(" << i+1 << "): "
          << info.error_info(i) << std::endl << std::endl;
  }
  return 0;
}
