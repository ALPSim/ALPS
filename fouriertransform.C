/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
*                       Matthias Troyer <troyer@comp-phys.org>
*                       Emanuel Gull <gullc@comp-phys.org>
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
/*fouriertransform.C*/
#include "fouriertransform.h"
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <alps/parameters.h>

void generate_spline_matrix(dense_matrix & spline_matrix) {

  // spline_matrix has dimension (N+1)x(N+1)
  int Np1 = spline_matrix.size1();

  // A is the matrix whose inverse defines spline_matrix
  //   
  //      6                   6
  //      1  4  1
  //         1  4  1
  // A =        ...
  //
  //                    1  4  1
  //     -2 -1             1  2    
  
  dense_matrix A = 4*boost::numeric::ublas::identity_matrix<double>(Np1);
  
  for (int i=1; i<Np1-1; i++) {
    A(i,i-1) = 1;
    A(i,i+1) = 1;
  }
  A(0,0) = 6;
  A(0, Np1-1) = 6;
  A(Np1-1, 0) = -2;
  A(Np1-1, 1) = -1;
  A(Np1-1, Np1-2) = 1;
  A(Np1-1, Np1-1) = 2;
  //std::cout<<A<<std::endl;  
  // solve A*spline_matrix=I
  // gesv solves A*X=B, input for B is I, output (=solution X) is spline_matrix
  /*spline_matrix = boost::numeric::ublas::identity_matrix<double>(Np1);   
  boost::numeric::bindings::lapack::gesv(A, spline_matrix);
  std::cout<<prec_prod(A, spline_matrix)<<std::endl;
*/
  //Emanuel avoid inversion of matrix
  spline_matrix=A;
}
void solve_second_derivatives(dense_matrix &spline_matrix, std::vector<double> &rhs){
  //avoid matrix inversion and instead solve the system of equations: S *v2 =
  //rhs for v2. The result will be found in rhs as we use the LAPACK call dgesv
  dense_matrix v(rhs.size(), 1);
  //std::cout<<v<<std::endl;
  for(uint i=0;i<rhs.size();++i){
    v(i,0)=rhs[i];
  }
  boost::numeric::bindings::lapack::gesv(spline_matrix, v);
  //std::cout<<rhs<<std::endl;
  for(uint i=0;i<rhs.size();++i){
    rhs[i]=v(i,0);
  }
}

///same as the template function in fouriertransform.h, but not for the
//multiple vectors. Use green_function_t instead.
void FourierTransformer::backward_ft(itime_green_function_t &G_tau, const matsubara_green_function_t &G_omega) const {
  assert(G_tau.nflavor()==G_omega.nflavor() && G_tau.nsite()==G_omega.nsite() && G_tau.nsite()==1);
  uint N_tau = G_tau.ntime()-1;
  uint N_omega = G_omega.nfreq();
  matsubara_green_function_t G_omega_no_model(G_omega);
  double dt = beta_/N_tau;
             
  for(uint f=0;f<G_omega.nflavor();++f){
    for (uint k=0; k<N_omega; k++) {
      std::complex<double> iw(0,(2*k+1)*M_PI/beta_);
      G_omega_no_model(k,f) -= f_omega(iw, c1_,c2_, c3_);
    }     
    for (uint i=0; i<N_tau; i++) {
      G_tau(i,f) = f_tau(i*dt, beta_, c1_, c2_, c3_); 
      for (uint k=0; k<N_omega; k++) {
        double wt((2*k+1)*i*M_PI/N_tau);
        G_tau(i,f) += 2/beta_*(cos(wt)*G_omega_no_model(k,f).real()+sin(wt)*G_omega_no_model(k,f).imag());
      }
    }
    G_tau(N_tau,f)=-1-G_tau(0,f);
  }
}
void FourierTransformer::forward_ft(const itime_green_function_t & G_tau, matsubara_green_function_t & G_omega) const {

  uint N=G_tau.ntime()-1;
  double dt = beta_/N;
  dense_matrix spline_matrix(N+1, N+1);
  assert(G_omega.nfreq()==N); //we implicitely assume that we have as many imag time slices as matsubara frequencies. If this is not the case: interpolate,
  std::vector<double> G_tau_vector(G_tau.ntime());
  std::vector<double> G_tau_derivatives(G_tau.ntime());
  assert(c1_==1); //otherwise change evaluate_second_derivatives!
  for(uint flavor=0;flavor<G_tau.nflavor();++flavor){
    generate_spline_matrix(spline_matrix);
    for(uint i=0;i<G_tau.ntime();++i){
      G_tau_vector[i]=G_tau(i,flavor);
    }
    evaluate_second_derivatives(dt, c2_, c3_, spline_matrix, G_tau_vector, G_tau_derivatives);
    for (uint k=0; k<N; k++) {
      std::complex<double> iw(0, (2*k+1)*M_PI/beta_);
      G_omega(k,flavor) = 0;
                                            
      for (uint n=1; n<N; n++) {
        G_omega(k,flavor) += exp(iw*(n*dt))*(G_tau_derivatives[n+1]-2*G_tau_derivatives[n]+G_tau_derivatives[n-1]);
      }
      G_omega(k,flavor) += G_tau_derivatives[1] - G_tau_derivatives[0] + G_tau_derivatives[N] - G_tau_derivatives[N-1];
      G_omega(k,flavor) *= 1./(dt*pow(iw,4));
      G_omega(k,flavor) += f_omega(iw, c1_, c2_, c3_);
    }
  }
} 
void FourierTransformer::generate_transformer(const alps::Parameters &parms,boost::shared_ptr<FourierTransformer> & fourier_ptr){
  if(parms.defined("GENERAL_FOURIER")){
    std::cout<<"using General Fourier transform"<<std::endl;
    fourier_ptr.reset(new FourierTransformer(parms["BETA"], parms["MU"]));
  } else{
    std::cout<<"using Bethe lattice Fourier transform"<<std::endl;
    fourier_ptr.reset(new BetheFourierTransformer(parms["BETA"], parms["MU"], parms["t"]));
  }
}

