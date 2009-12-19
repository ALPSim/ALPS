/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2004 by Salvatore R. Manmana <Salva@theo3.physik.uni-stuttgart.de>,
*                            Reinhard M. Noack <Reinhard.Noack@physik.uni-marburg.de>
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

#include "densmatrix.h"

DensityMatrix::DensityMatrix(const WaveFunction& psi, LR lr) 
{

  int imax = psi.size1();
  int jmax = psi.size2();  
  
  if (lr == Left)
    {
      rho.resize(imax,imax);
      rho = boost::numeric::ublas::prec_prod(psi,trans(psi));
    }
  else
    {
      rho.resize(jmax,jmax);
      rho = boost::numeric::ublas::prec_prod(trans(psi),psi);
    }      
}

void DensityMatrix::newBasis(int cutoff, Matrix& basis, double& discarded_weight)
{
  using namespace boost::numeric::ublas;
  Vector singvals(rho.size1());
  
  //UBLAS BINDINGS! -> which header-file??
  //boost::numeric::bindings::lapack::syev('V','U',rho,singvals);
  ietl2lapack::syev(rho.size1(),&(rho(0,0)),&(singvals(0)),'V','U');
  
  if(cutoff < rho.size1()){
    
    basis.resize(rho.size1(),cutoff);
    basis.clear();

    // NOTE: the LAPACK-routine dsyev returns the eigenvalues of the matrix
    // in ascending order; the eigenvectors are ordered correspondingly. 
    // This means here, that the LAST eigenvectors have the largest
    // density-matrix eigenvalues!
    
    // rewrite using iterators:
    for(int j = 0; j < cutoff; j++){
      int ev_numb = rho.size2() - 1 - j;
      matrix_column<Matrix>(basis,j) = matrix_column<Matrix>(rho,ev_numb);
    }
    
    // calculate discarded weight:
    vector_range<Vector> kept_singvals (singvals, range (singvals.size() - cutoff, singvals.size()));
    discarded_weight = 1.0 - sum(kept_singvals);
    
    cout << "sum over all singular values = " << sum(singvals) <<endl;
    cout << "discarded weight = " << discarded_weight <<endl;
  }
  else{
    
    basis = rho;
    discarded_weight = 1.0 - sum(singvals);
    
    cout 
      << "keeping all density matrix states; calculated discarded weight = " 
      << discarded_weight <<endl;
  
  }
}
