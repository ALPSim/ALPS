/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
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


#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>

#include "fouriertransform.h"
#include <valarray>
#include <boost/numeric/ublas/io.hpp>
#include <alps/parameter.h>

#include <math.h>

typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;


void FourierTransformer::backward_ft(itime_green_function_t &G_tau, 
                                     const matsubara_green_function_t &G_omega) const {
  assert(G_tau.nflavor()==G_omega.nflavor() && G_tau.nsite()==G_omega.nsite());
  unsigned int N_tau = G_tau.ntime()-1;
  unsigned int N_omega = G_omega.nfreq();
  unsigned int N_site = G_omega.nsite();
  matsubara_green_function_t G_omega_no_model(G_omega);
  matsubara_green_function_t test(G_omega);
  double dt = beta_/N_tau;
  for(unsigned int f=0;f<G_omega.nflavor();++f){
    for (unsigned int s1=0; s1<N_site; ++s1){
      for (unsigned int s2=0; s2<N_site; ++s2) {
        if(c1_[f][s1][s2]==0 && c2_[f][s1][s2] == 0 && c3_[f][s1][s2]){  //nothing happening in this gf.
          for (unsigned int i=0; i<=N_tau; i++) {
            G_tau(i,s1,s2,f)=0.;
          }
        } 
        else {
          for (unsigned int k=0; k<N_omega; k++) {
            std::complex<double> iw(0,(2*k+1)*M_PI/beta_);
            G_omega_no_model(k,s1,s2,f) -= f_omega(iw, c1_[f][s1][s2],c2_[f][s1][s2], c3_[f][s1][s2]);
          } 
          for (unsigned int i=0; i<N_tau; i++) {
            G_tau(i,s1,s2,f) = f_tau(i*dt, beta_, c1_[f][s1][s2], c2_[f][s1][s2], c3_[f][s1][s2]); 
            for (unsigned int k=0; k<N_omega; k++) {
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
}


void generate_spline_matrix(dense_matrix & spline_matrix, double dt) {
  
  // spline_matrix has dimension (N+1)x(N+1)
  int Np1 = spline_matrix.size1();
  //std::cout<<"spline matrix size is: "<<Np1<<std::endl;
  // A is the matrix whose inverse defines spline_matrix
  //   
  //      6                   6
  //      1  4  1
  //         1  4  1
  // A =        ...
  //
  //                    1  4  1
  //     -2 -1     0       1  2
  spline_matrix.clear(); 
  dense_matrix A = 4*dt/6.*boost::numeric::ublas::identity_matrix<double>(Np1);
  
  for (int i=1; i<Np1-1; i++) {
    A(i,i-1) = dt/6.;
    A(i,i+1) = dt/6.;
  }
  A(0,0) = 1.;
  A(0, Np1-1) = 1.;
  A(Np1-1, 0) = -2.*dt/6.;
  A(Np1-1, 1) = -1.*dt/6.;
  A(Np1-1, Np1-2) = 1*dt/6.;
  A(Np1-1, Np1-1) = 2*dt/6.;
  
  // solve A*spline_matrix=I
  // gesv solves A*X=B, input for B is I, output (=solution X) is spline_matrix
  spline_matrix = boost::numeric::ublas::identity_matrix<double>(Np1);   
  boost::numeric::ublas::vector<fortran_int_t> ipivot(A.size1());
  boost::numeric::bindings::lapack::gesv(A, ipivot,spline_matrix);
}



void evaluate_second_derivatives(double dt/*, double BETA*/, dense_matrix & spline_matrix, std::vector<double> & g, std::vector<double> & second_derivatives, const double c1g, const double c2g, const double c3g) {
  
  // g, rhs and second_derivatives have dimension N+1
  int Np1 = spline_matrix.size1();
  //assert(c1g==1); 
  // rhs is the vector containing the data of the curve y = g(tau), which allows to 
  // compute the vector of second derivatives y'' at times tau_n by evaluating
  // y'' = spline_matrix * rhs(y)
  //
  //                         0                                
  //                         y0 - 2*y1 + y2
  //                         y1 - 2*y2 + y3
  // rhs = 6/(delta_tau)^2 * ...
  //
  //                         yNp1-3 - 2*yNp1-2 + yNp1-1
  //                         y0 - y1 + yNp1-2 - yNp1-1     
  
  std::vector<double> rhs(Np1, 0);
  std::cout<<"constants: "<<c1g<<" "<<c2g<<" "<<c3g<<std::endl;
  rhs[0] = -c3g; //G''(0)+G''(beta)=-c3
  for (int i=1; i<Np1-1; i++) {
    rhs[i] = (g[i-1]-2*g[i]+g[i+1])/dt;
  }
  rhs[Np1-1] = c2g -1./dt*(-g[0] + g[1] -g[Np1-2] + g[Np1-1]);
  
  for (int i=0; i<Np1; i++) {
    second_derivatives[i]=0;
    for (int j=0; j<Np1; j++) {
      second_derivatives[i] += spline_matrix(i,j)*rhs[j];
    }
  }
}



void FourierTransformer::forward_ft(const itime_green_function_t & gtau, matsubara_green_function_t & gomega) const 
{
  std::vector<double> v(gtau.ntime());
  std::vector<std::complex<double> > v_omega(gomega.nfreq());
  int Np1 = v.size();
  int N = Np1-1;
  int N_omega = v_omega.size();
  double dt = beta_/N;
  
  for(unsigned f=0;f<gtau.nflavor();++f){
    for(unsigned p=0;p<gtau.nsite();++p){
      for(unsigned q=0;q<gtau.nsite();++q){
        for(int tau=0;tau<Np1;++tau){
          v[tau]=gtau(tau,p,q,f);
        }
        
        dense_matrix spline_matrix(Np1, Np1);
        generate_spline_matrix(spline_matrix, dt);
        // matrix containing the second derivatives y'' of interpolated y=v[tau] at points tau_n 
        std::vector<double> v2(Np1, 0); 
        evaluate_second_derivatives(dt/*,beta_*/, spline_matrix, v, v2, c1_[f][p][q], c2_[f][p][q], c3_[f][p][q]);
        v_omega.assign(N_omega, 0);
        
        for (int k=0; k<N_omega; k++) {
          std::complex<double> iw(0, M_PI*(2*k+1)/beta_);
          for (int n=1; n<N; n++) {
            v_omega[k] += exp(iw*(n*dt))*(v2[n+1]-2*v2[n]+v2[n-1]); //partial integration, four times. Then approximate the fourth derivative by finite differences
          }
          v_omega[k] += (v2[1] - v2[0] + v2[N] - v2[N-1]); //that's the third derivative, on the boundary
          v_omega[k] *= 1./(dt*iw*iw*iw*iw);
          v_omega[k] += f_omega(iw, c1_[f][p][q], c2_[f][p][q], c3_[f][p][q]); //that's the boundary terms of the first, second, and third partial integration.
          //std::cout<<"derivative at boundary: "<<-v2[1] + v2[0] + v2[N] - v2[N-1]<<" divisor: "<< 1./(dt*iw*iw*iw*iw)<<std::endl;
          gomega(k,p,q,f)=v_omega[k]; //that's the proper convention for the self consistency loop here.
        }
      }
    }
  }
} 



void FourierTransformer::append_tail(matsubara_green_function_t& G_omega, 
                                     const matsubara_green_function_t& G0_omega,
                                     const int nfreq_measured) const
{
  for (spin_t flavor=0; flavor<G0_omega.nflavor(); ++flavor) {
    for (site_t k=0; k<G0_omega.nsite(); ++k) {
      std::cout << "append tail to self-energy with coefficients: "
      << " " << Sc0_[flavor][k][k]
      << " " << Sc1_[flavor][k][k]
      << " " << Sc2_[flavor][k][k] << std::endl;
      for (frequency_t freq=nfreq_measured; freq<G0_omega.nfreq(); ++freq) {
        std::complex<double> iw(0,(2*freq+1)*M_PI/beta_);
        std::complex<double> Sigma = Sc0_[flavor][k][k] + Sc1_[flavor][k][k]/iw + Sc2_[flavor][k][k]/(iw*iw);
        G_omega(freq, k, k, flavor) = 1./(1./G0_omega(freq, k, k, flavor) - Sigma);
      }
    }
  }
}

void FourierTransformer::generate_transformer(const alps::Parameters &parms,
                                              boost::shared_ptr<FourierTransformer> &fourier_ptr)
{
  int n_flavors = parms.value_or_default("FLAVORS", 2); 
  if (static_cast<int>(parms.value_or_default("SITES", 1))!=1) 
    throw std::logic_error("ERROR: FourierTransformer::generate_transformer : SITES!=1, for cluster fourier transforms please use the cluster version of this framework");
  double h = static_cast<double>(parms.value_or_default("H",0.));
  std::cout << "using general fourier transformer" << "\n";
  std::vector<double> eps(n_flavors);
  std::vector<double> epssq(n_flavors);
  for (int f=0; f<n_flavors; ++f) {
    eps[f] = parms.value_or_default("EPS_"+boost::lexical_cast<std::string>(f),0.0);
    epssq[f] = parms.value_or_default("EPSSQ_"+boost::lexical_cast<std::string>(f),1.0);
  }
  fourier_ptr.reset(new SimpleG0FourierTransformer((double)parms["BETA"], (double)parms["MU"], h,
                                                   n_flavors, eps, epssq));
}




void FourierTransformer::generate_transformer_U(const alps::Parameters &parms,
                                                boost::shared_ptr<FourierTransformer> &fourier_ptr,
                                                const std::vector<double> &densities)
{
  int n_flavors = parms.value_or_default("FLAVORS", 2); 
  if (static_cast<int>(parms.value_or_default("SITES", 1))!=1) 
    throw std::logic_error("ERROR: FourierTransformer::generate_transformer_U : SITES!=1, for cluster fourier transforms please use the cluster version of this framework");
  double U = parms["U"];
  std::cout << "using general fourier transformer" << "\n";
  std::vector<std::vector<double> > eps(n_flavors,std::vector<double>(1));
  std::vector<std::vector<double> >epssq(n_flavors,std::vector<double>(1));
  for (int f=0; f<n_flavors; ++f) {
    eps[f][0] = parms.value_or_default("EPS_"+boost::lexical_cast<std::string>(f),0.0);
    epssq[f][0] = parms.value_or_default("EPSSQ_"+boost::lexical_cast<std::string>(f),1.0);
  }
  fourier_ptr.reset(new GFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., U, 
                                            n_flavors, 1, densities, eps, epssq));
}
