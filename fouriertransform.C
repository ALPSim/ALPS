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
#include <valarray>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <alps/parameter.h>

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
void FourierTransformer::backward_ft(itime_green_function_t &G_tau, 
				     const matsubara_green_function_t &G_omega) const {
  assert(G_tau.nflavor()==G_omega.nflavor() && G_tau.nsite()==G_omega.nsite());
  uint N_tau = G_tau.ntime()-1;
  uint N_omega = G_omega.nfreq();
  uint N_site = G_omega.nsite();
  matsubara_green_function_t G_omega_no_model(G_omega);
  double dt = beta_/N_tau;
  for(uint f=0;f<G_omega.nflavor();++f){
    for (uint s1=0; s1<N_site; ++s1){
      for (uint s2=0; s2<N_site; ++s2) {
	if (s1==s2) {
	  for (uint k=0; k<N_omega; k++) {
	    std::complex<double> iw(0,(2*k+1)*M_PI/beta_);
	    //std::cout<<G_omega_no_model(k,f)<<std::endl;
	    G_omega_no_model(k,s1,s2,f) -= f_omega(iw, c1_,c2_, c3_);
	    //std::cout<<k<<" "<<G_omega_no_model(k,f)<<std::endl;
	  } 
	}
	//std::cout<<c1_<<" "<<c2_<<" "<<c3_<<std::endl;
	for (uint i=0; i<N_tau; i++) {
	  G_tau(i,s1,s2,f) = s1==s2 ? f_tau(i*dt, beta_, c1_, c2_, c3_) : 0.; 
	  //std::cout<<i<<" "<<G_tau(i,f)<<std::endl;
	  for (uint k=0; k<N_omega; k++) {
	    double wt((2*k+1)*i*M_PI/N_tau);
	    G_tau(i,s1,s2,f) += 2/beta_*(cos(wt)*G_omega_no_model(k,s1,s2,f).real()+
					 sin(wt)*G_omega_no_model(k,s1,s2,f).imag());
	  }
	}
	G_tau(N_tau,s1,s2,f)= s1==s2 ? -1. : 0.;
	G_tau(N_tau,s1,s2,f)-=G_tau(0,s1,s2,f);
      }
    }
  }
}


void FourierTransformer::forward_ft(const itime_green_function_t & G_tau, matsubara_green_function_t & G_omega) const {
  assert(G_tau.nflavor()==G_omega.nflavor() && G_tau.nsite()==G_omega.nsite());
  uint N=G_tau.ntime()-1;
  uint N_site=G_tau.nsite();
  //double dt = beta_/N;
  uint n_oversampling_t=10; //not sure that helps... check again!
  uint n_oversampling_w=1; //not sure that helps... check again!
  uint N_large_t=(G_tau.ntime()-1)*n_oversampling_t+1;
  uint N_large_w=(G_tau.ntime()-1)*n_oversampling_w+1;
  double dt=beta_/(N_large_t-1);
//  assert(G_omega.nfreq()==N); 
  //we implicitely assume that we have as many imag time slices as matsubara frequencies. 
  //If this is not the case: interpolate,
  double G_tau_vector[N_large_t];
  std::valarray<std::complex<double> > G_omega_vector(N_large_w);
  //std::vector<double> G_tau_derivatives(G_tau.ntime());
  assert(c1_==1); //otherwise change evaluate_second_derivatives!
  for(uint flavor=0;flavor<G_tau.nflavor();++flavor){
    for (uint s1=0; s1<N_site; ++s1){
      for (uint s2=0; s2<N_site; ++s2) {
	for(uint i=0;i<N;++i){
	  for(uint j=0;j<n_oversampling_t;++j){
	    G_tau_vector[i*n_oversampling_t+j]=
	      (1./n_oversampling_t*(G_tau(i,s1,s2,flavor)*(n_oversampling_t-j)+G_tau(i+1,s1,s2,flavor)*j));
	  }
	  G_tau_vector[N_large_t-1]=G_tau(N,s1,s2,flavor);
	}
	for (uint n=0; n<N_large_w; n++) {
	  double w_n=(2*n+1)*M_PI/beta_;
	  G_omega_vector[n]=0;
	  for(uint i=0;i<N_large_t-1;++i){
	    //if(n==0)
	    //std::cout<<i<<" "<<dt<<" "<<i*dt<<std::endl;
            double phi=w_n*i*dt;
	    if (s1==s2) 
	      G_omega_vector[n] += std::complex<double>(cos(phi), sin(phi))*
		((G_tau_vector[i]-f_tau(dt*i,beta_, c1_, c2_, c3_))*dt);
	    else
	      G_omega_vector[n] += dt*std::exp((std::complex<double>(0,1))*(w_n*i*dt))*G_tau_vector[i];
	  }
	  //std::cerr<<n<<" "<<G_omega_vector[n].real()<<" "<<G_omega_vector[n].imag()<<std::endl;
	}
        std::cout<<"doing forward ft for nfreq= "<<G_omega.nfreq()<<std::endl;
	for(uint n=0;n<G_omega.nfreq();++n){
	  G_omega(n,s1,s2,flavor)=n<N_large_w?G_omega_vector[n]:0.;
	  if (s1==s2) {
	    std::complex<double> iw(0,(2*n+1)*M_PI/beta_);
	    G_omega(n,s1,s2,flavor) += f_omega(iw, c1_,c2_, c3_);
	  }
	}
      }
    }
  }
}

 
void FourierTransformer::generate_transformer(const alps::Parameters &parms,boost::shared_ptr<FourierTransformer> & fourier_ptr, double mu_shift){
  if(parms.defined("GENERAL_FOURIER") || parms.defined("CLUSTER_LOOP")){
    std::cout<<"using General Fourier transform"<<std::endl;
    fourier_ptr.reset(new FourierTransformer(parms["BETA"], (double)(parms["MU"])+mu_shift));
  } else{
    std::cout<<"using Bethe lattice Fourier transform"<<std::endl;
    fourier_ptr.reset(new BetheFourierTransformer(parms["BETA"], (double)(parms["MU"])+mu_shift, parms["t"]));
  }
}


void FourierTransformer::generate_transformer_U(const alps::Parameters &parms,boost::shared_ptr<FourierTransformer> & fourier_ptr){
  if(parms.defined("GENERAL_FOURIER") || parms.defined("CLUSTER_LOOP")){
    std::cout<<"using General Fourier transform"<<std::endl;
    fourier_ptr.reset(new FourierTransformer(parms["BETA"], parms["MU"]));
  } else{
    std::cout<<"using Bethe lattice Fourier transform"<<std::endl;
    fourier_ptr.reset(new BetheFourierTransformer(parms["BETA"], parms["MU"], parms["t"], parms["U"]));
  }
}
