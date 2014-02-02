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
#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>

#include "fouriertransform.h"
#include "cluster.h"
#include <valarray>
#include <alps/parameter.h>



///same as the template function in fouriertransform.h, but not for the
//multiple vectors. Use green_function_t instead.
void FourierTransformer::frequency_to_time_ft(itime_green_function_t &G_tau,
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
        if(G_omega.shape()==diagonal && (s1 != s2)) continue;
        if(G_omega.shape()==blockdiagonal){
          bool is_nonzero=false;
          if(s1==s2) is_nonzero=true;
          if(s1%2==0 && s2==s1+1) is_nonzero=true;
          if(s1%2==1 && s2==s1-1) is_nonzero=true;
          if (!is_nonzero){  //nothing happening in this gf.
            /*for (unsigned int i=0; i<=N_tau; i++) {
             G_tau(i,s1,s2,f)=0.;
             }*/
            continue;
          }
        }
        //this is the standard case: we need to do a FT.
        {
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

void FourierTransformer::frequency_to_time_ft(nambu_itime_green_function_t &G_tau,
                                              const nambu_matsubara_green_function_t &G_omega) const {
  assert(G_tau.nflavor()==G_omega.nflavor() && G_tau.nsite()==G_omega.nsite());
  unsigned int N_tau = G_tau.ntime()-1;
  unsigned int N_omega = G_omega.nfreq();
  unsigned int N_site = G_omega.nsite();
  nambu_matsubara_green_function_t G_omega_no_model(G_omega);
  double dt = beta_/N_tau;
  std::vector<std::vector<double> >cos_wt(N_tau, std::vector<double>(N_omega));
  std::vector<std::vector<double> >sin_wt(N_tau, std::vector<double>(N_omega));
  for (unsigned int i=0; i<N_tau; i++) {
    for (unsigned int k=0; k<N_omega; k++) {
      double wt((2*k+1)*i*M_PI/N_tau);
      cos_wt[i][k]=cos(wt);
      sin_wt[i][k]=sin(wt);
    }
  }
  for (unsigned int s1=0; s1<N_site; ++s1){
    for (unsigned int s2=0; s2<N_site; ++s2) {
      if((c1_[0][s1][s2]==0) && (c2_[0][s1][s2] == 0) && (c3_[0][s1][s2]== 0)){  //nothing happening in this gf - this fct is completely zero.
        for (unsigned int i=0; i<=N_tau; i++) {
          G_tau(i,s1,s2,0,0)=0.;
          G_tau(i,s1,s2,0,1)=0.;
          G_tau(i,s1,s2,1,0)=0.;
          G_tau(i,s1,s2,1,1)=0.;
        }
      }
      else {
        for (unsigned int k=0; k<N_omega; k++) {
          std::complex<double> iw(0,(2*k+1)*M_PI/beta_);
          G_omega_no_model(k,s1,s2,0,0) -= f_omega(iw, c1_[0][s1][s2],c2_[0][s1][s2], c3_[0][s1][s2]);
          G_omega_no_model(k,s1,s2,0,1) -= 0;
          G_omega_no_model(k,s1,s2,1,0) -= 0;
          G_omega_no_model(k,s1,s2,1,1) -= f_omega(iw, c1_[1][s1][s2],-c2_[1][s1][s2], c3_[1][s1][s2]); //the '-' here is because we have Nambu notation.
          
          //std::cout<<k<<" "<<G_omega_no_model(k,s1,s2,0,0).real()<<" "<<G_omega_no_model(k,s1,s2,0,0).imag()<<" "<<G_omega_no_model(k,s1,s2,1,1).real()<<" "<<G_omega_no_model(k,s1,s2,1,1).imag()<<std::endl;
        }
        for (unsigned int i=0; i<N_tau; i++) {
          G_tau(i,s1,s2,0,0) = f_tau(i*dt, beta_, c1_[0][s1][s2], c2_[0][s1][s2], c3_[0][s1][s2]);
          G_tau(i,s1,s2,0,1) = 0;
          G_tau(i,s1,s2,1,0) = 0;
          G_tau(i,s1,s2,1,1) = f_tau(i*dt, beta_, c1_[1][s1][s2], -c2_[1][s1][s2], c3_[1][s1][s2]);
          for (unsigned int k=0; k<N_omega; k++) {
            G_tau(i,s1,s2,0,0) += 2/beta_*(cos_wt[i][k]*G_omega_no_model(k,s1,s2,0,0).real()+sin_wt[i][k]*G_omega_no_model(k,s1,s2,0,0).imag());
            G_tau(i,s1,s2,0,1) += 2/beta_*(cos_wt[i][k]*G_omega_no_model(k,s1,s2,0,1).real()+sin_wt[i][k]*G_omega_no_model(k,s1,s2,0,1).imag());
            G_tau(i,s1,s2,1,0) += 2/beta_*(cos_wt[i][k]*G_omega_no_model(k,s1,s2,1,0).real()+sin_wt[i][k]*G_omega_no_model(k,s1,s2,1,0).imag());
            G_tau(i,s1,s2,1,1) += 2/beta_*(cos_wt[i][k]*G_omega_no_model(k,s1,s2,1,1).real()+sin_wt[i][k]*G_omega_no_model(k,s1,s2,1,1).imag());
          }
        }
        G_tau(N_tau,s1,s2,0,0)= (s1==s2 ? -1. : 0.) - G_tau(0,s1,s2,0,0);
        G_tau(N_tau,s1,s2,0,1)=  - G_tau(0,s1,s2,0,1);
        G_tau(N_tau,s1,s2,1,0)=  - G_tau(0,s1,s2,1,0);
        G_tau(N_tau,s1,s2,1,1)= (s1==s2 ? -1. : 0.) - G_tau(0,s1,s2,1,1);
      }
      for(unsigned int i=0;i<=N_tau;++i){
        //panic... this means we have a real problem.
        if(G_tau(i,s1,s2,0,0)<=-1.01 || G_tau(i,s1,s2,0,0) >=1.01){ std::cout<<G_tau(i,s1,s2,0,0)<<" (s1,s2,f1,f2): "<<s1<<" "<<s2<<" "<<0<<" "<<0<<std::endl; throw std::runtime_error("problem with G0 at index."); }
        if(G_tau(i,s1,s2,0,1)<=-1.01 || G_tau(i,s1,s2,0,1) >=1.01){ std::cout<<G_tau(i,s1,s2,0,1)<<" (s1,s2,f1,f2): "<<s1<<" "<<s2<<" "<<0<<" "<<1<<std::endl; throw std::runtime_error("problem with G0 at index."); }
        if(G_tau(i,s1,s2,1,0)<=-1.01 || G_tau(i,s1,s2,1,0) >=1.01){ std::cout<<G_tau(i,s1,s2,1,0)<<" (s1,s2,f1,f2): "<<s1<<" "<<s2<<" "<<1<<" "<<0<<std::endl; throw std::runtime_error("problem with G0 at index."); }
        if(G_tau(i,s1,s2,1,1)<=-1.01 || G_tau(i,s1,s2,1,1) >=1.01){ std::cout<<G_tau(i,s1,s2,1,1)<<" (s1,s2,f1,f2): "<<s1<<" "<<s2<<" "<<1<<" "<<1<<std::endl; throw std::runtime_error("problem with G0 at index."); }
      }
      for(unsigned int i=0;i<=N_tau;++i){
        //roundoff error. We should probably use more frequencies.
        if(G_tau(i,s1,s2,0,0)<=-1.0){ G_tau(i,s1,s2,0,0)=-1.0+1.e-6; std::cerr<<"adjusting G at "<<i<<" "<<s1<<" "<<s2<<" "<<0<<" "<<0<<" to be larger than -1"<<std::endl; }
        if(G_tau(i,s1,s2,0,1)<=-1.0){ G_tau(i,s1,s2,0,1)=-1.0+1.e-6; std::cerr<<"adjusting G at "<<i<<" "<<s1<<" "<<s2<<" "<<0<<" "<<1<<" to be larger than -1"<<std::endl; }
        if(G_tau(i,s1,s2,1,0)<=-1.0){ G_tau(i,s1,s2,1,0)=-1.0+1.e-6; std::cerr<<"adjusting G at "<<i<<" "<<s1<<" "<<s2<<" "<<1<<" "<<0<<" to be larger than -1"<<std::endl; }
        if(G_tau(i,s1,s2,1,1)<=-1.0){ G_tau(i,s1,s2,1,1)=-1.0+1.e-6; std::cerr<<"adjusting G at "<<i<<" "<<s1<<" "<<s2<<" "<<1<<" "<<1<<" to be larger than -1"<<std::endl; }
      }
    }
  }
}





void generate_spline_matrix(dense_matrix & spline_matrix, double dt) {
  
  // spline_matrix has dimension (N+1)x(N+1)
  int Np1 = spline_matrix.size1();
  std::cout<<"spline matrix size is: "<<Np1<<std::endl;
  // A is the matrix whose inverse defines spline_matrix
  //
  //      6                   6
  //      1  4  1
  //         1  4  1
  // A =        ...
  //
  //                    1  4  1
  //    -2  0              0  2
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



void FourierTransformer::time_to_frequency_ft(const itime_green_function_t & gtau, matsubara_green_function_t & gomega) const
{
  std::vector<double> v(gtau.ntime());
  std::vector<std::complex<double> > v_omega(gomega.nfreq());
  int Np1 = v.size();
  int N = Np1-1;
  int N_omega = v_omega.size();
  double dt = beta_/N;
  
  for(unsigned int f=0;f<gtau.nflavor();++f){
    for(unsigned int p=0;p<gtau.nsite();++p){
      for(unsigned int q=0;q<gtau.nsite();++q){
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
          gomega(k,p,q,f)=v_omega[k]; //that's the proper convention for the self consistency loop here.
        }
        
      }
    }
  }
}

void FourierTransformer::time_to_frequency_ft(const nambu_itime_green_function_t & gtau, nambu_matsubara_green_function_t & gomega) const
{
  std::vector<double> v(gtau.ntime());
  std::vector<std::complex<double> > v_omega(gomega.nfreq());
  int Np1 = v.size();
  int N = Np1-1;
  int N_omega = v_omega.size();
  double dt = beta_/N;
  
  for(unsigned int f1=0;f1<gtau.nflavor();++f1){
    for(unsigned int f2=0;f2<gtau.nflavor();++f2){
      for(unsigned int p=0;p<gtau.nsite();++p){
        for(unsigned int q=0;q<gtau.nsite();++q){
          for(int tau=0;tau<Np1;++tau){
            v[tau]=gtau(tau,p,q,f1,f2);
          }
          
          dense_matrix spline_matrix(Np1, Np1);
          generate_spline_matrix(spline_matrix, dt);
          // matrix containing the second derivatives y'' of interpolated y=v[tau] at points tau_n
          std::vector<double> v2(Np1, 0);
          if(f1==f2){
            evaluate_second_derivatives(dt, spline_matrix, v, v2, c1_[f1][p][q], c2_[f1][p][q], c3_[f1][p][q]);
          }else{
            evaluate_second_derivatives(dt, spline_matrix, v, v2, 0, 0, 0);
          }
          v_omega.assign(N_omega, 0);
          
          for (int k=0; k<N_omega; k++) {
            std::complex<double> iw(0, M_PI*(2*k+1)/beta_);
            for (int n=1; n<N; n++) {
              v_omega[k] += exp(iw*(n*dt))*(v2[n+1]-2*v2[n]+v2[n-1]); //partial integration, four times. Then approximate the fourth derivative by finite differences
            }
            v_omega[k] += (v2[1] - v2[0] + v2[N] - v2[N-1]); //that's the third derivative, on the boundary
            v_omega[k] *= 1./(dt*iw*iw*iw*iw);
            if(f1==f2){
              v_omega[k] += f_omega(iw, c1_[f1][p][q], c2_[f1][p][q], c3_[f1][p][q]); //that's the boundary terms of the first, second, and third partial integration.
            }
            gomega(k,p,q,f1,f2)=v_omega[k]; //that's the proper convention for the self consistency loop here.
          }
        }
      }
    }
  }
}


void FourierTransformer::append_tail(matsubara_green_function_t& G_omega,
                                     const matsubara_green_function_t& G0_omega,
                                     const unsigned int nfreq_measured) const
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



void ClusterAntiFerroGFourierTransformer::append_tail(matsubara_green_function_t& G_omega,
                                                      const matsubara_green_function_t& G0_omega,
                                                      const unsigned int nfreq_measured) const
{
  matsubara_green_function_t selfenergy_tail(G0_omega.nfreq()-nfreq_measured, G0_omega.nsite(), G0_omega.nflavor(), G_omega.shape());
  selfenergy_tail.clear();
  matsubara_green_function_t G_omega_tail(selfenergy_tail);
  for (spin_t flavor=0; flavor<G0_omega.nflavor(); ++flavor) {
    for (frequency_t freq=nfreq_measured; freq<G0_omega.nfreq(); ++freq) {
      std::complex<double> iw(0,(2*freq+1)*M_PI/beta_);
      for (site_t k=0; k<G0_omega.nsite(); k+=2) {
        for (int p1=0; p1<=1; p1++) {
          for (int p2=0; p2<=1; p2++) {
            selfenergy_tail(freq-nfreq_measured, k+p1, k+p2, flavor) =
            Sc0_[flavor][k+p1][k+p2] + Sc1_[flavor][k+p1][k+p2]/iw + Sc2_[flavor][k+p1][k+p2]/(iw*iw);
            G_omega_tail(freq-nfreq_measured, k+p1, k+p2, flavor) = G0_omega(freq, k+p1, k+p2, flavor);
          }
        }
      }
    }
  }
  matsubara_green_function_t selfenergy_measured(nfreq_measured, G0_omega.nsite(), G0_omega.nflavor(), G_omega.shape());
  matsubara_green_function_t G_omega_measured(nfreq_measured, G0_omega.nsite(), G0_omega.nflavor(), G_omega.shape());
  selfenergy_measured.clear();
  for (spin_t flavor=0; flavor<G0_omega.nflavor(); ++flavor) {
    for (frequency_t freq=0; freq<nfreq_measured; ++freq) {
      for (site_t k=0; k<G0_omega.nsite(); k+=2) {
        for (int p1=0; p1<=1; p1++) {
          for (int p2=0; p2<=1; p2++) {
            selfenergy_measured(freq, k+p1, k+p2, flavor) = G0_omega(freq, k+p1, k+p2, flavor);
            G_omega_measured(freq, k+p1, k+p2, flavor) = G_omega(freq, k+p1, k+p2, flavor);
          }
        }
      }
    }
  }
  selfenergy_measured = AntiFerroDCATransformer::green_function_inverse(selfenergy_measured);
  selfenergy_measured =
  AntiFerroDCATransformer::green_function_subtract(selfenergy_measured,
                                                   AntiFerroDCATransformer::green_function_inverse(G_omega_measured));
  G_omega_tail = AntiFerroDCATransformer::green_function_inverse(G_omega_tail);
  G_omega_tail = AntiFerroDCATransformer::green_function_subtract(G_omega_tail, selfenergy_tail);
  G_omega_tail = AntiFerroDCATransformer::green_function_inverse(G_omega_tail);
  for (spin_t flavor=0; flavor<G0_omega.nflavor(); ++flavor){
    for (frequency_t freq=nfreq_measured; freq<G0_omega.nfreq(); ++freq){
      for (site_t k=0; k<G0_omega.nsite(); k+=2){
        for (int p1=0; p1<=1; p1++){
          for (int p2=0; p2<=1; p2++){
            G_omega(freq, k+p1, k+p2, flavor) = G_omega_tail(freq-nfreq_measured, k+p1, k+p2, flavor);
          }
        }
      }
    }
  }
}




void FourierTransformer::generate_transformer(const alps::Parameters &parms,
                                              boost::shared_ptr<FourierTransformer> &fourier_ptr)
{
  int n_flavors = parms.value_or_default("FLAVORS", 2);
  int n_site = parms.value_or_default("SITES", 1);
  if (parms.defined("GENERAL_FOURIER_TRANSFORMER")) {
    std::cout << "using general fourier transformer" << "\n";
    std::vector<double> eps(n_flavors);
    std::vector<double> epssq(n_flavors);
    for (int f=0; f<n_flavors; ++f) {
      eps[f] = parms["EPS_"+boost::lexical_cast<std::string>(f)];
      epssq[f] = parms["EPSSQ_"+boost::lexical_cast<std::string>(f)];
    }
    fourier_ptr.reset(new SimpleG0FourierTransformer((double)parms["BETA"], (double)parms["MU"], (double)parms["H"],
                                                     n_flavors, eps, epssq));
  }
  else if (parms.defined("CLUSTER_LOOP")){
    if(std::string(parms["CLUSTER_LOOP"])==std::string("DCA")){
      boost::shared_ptr<DCATransformer> clusterhandler;
      DCATransformer::generate_transformer(parms, clusterhandler);
      std::vector<std::vector<double> > eps_f(n_flavors, clusterhandler->epsav());
      std::vector<std::vector<double> > epssq_f(n_flavors, clusterhandler->epssqav());
      fourier_ptr.reset(new ClusterG0FourierTransformer((double)parms["BETA"], (double)parms["MU"], (double)parms["H"],
                                                        n_flavors, n_site, eps_f, epssq_f));
    }else if(std::string(parms["CLUSTER_LOOP"])==std::string("CDMFT")){
      std::cout<<"using CDMFT cluster Fourier transformer"<<std::endl;
      std::cout<<"warning: G0 moments are not implemented for CDMFT"<<std::endl;
      std::vector<std::vector<double> > eps_f(n_flavors, std::vector<double>(n_site*n_site,0.00));
      std::vector<std::vector<double> > epssq_f(n_flavors, std::vector<double>(n_site*n_site,0.0001));
      fourier_ptr.reset(new ClusterG0FourierTransformer((double)parms["BETA"], (double)parms["MU"], (double)parms["H"],
                                                        n_flavors, n_site, eps_f, epssq_f));
    }
    else throw std::logic_error("please figure out how to do fourier transforms in your cluster framwork if it is not DCA (2)");
  }
  else {
    std::cout << "using Bethe lattice fourier transform" << "\n";
    std::vector<double> eps(n_flavors);
    std::vector<double> epssq(n_flavors);
    double t = parms["t"];
    for (int f=0; f<n_flavors; ++f) {
      eps[f] = 0.;
      epssq[f] = t*t;
    }
    fourier_ptr.reset(new SimpleG0FourierTransformer((double)parms["BETA"], (double)parms["MU"], (double)parms["H"],
                                                     n_flavors, eps, epssq));
  }
}




void FourierTransformer::generate_transformer_U(const alps::Parameters &parms,
                                                boost::shared_ptr<FourierTransformer> &fourier_ptr,
                                                const std::vector<double> &densities)
{
  int n_flavors = parms.value_or_default("FLAVORS", 2);
  int n_site = parms.value_or_default("SITES", 1);
  double U = parms["U"];
  if (parms.defined("GENERAL_FOURIER_TRANSFORMER")) {
    std::cout << "using general fourier transformer" << "\n";
    std::vector<std::vector<double> > eps(n_flavors);
    std::vector<std::vector<double> >epssq(n_flavors);
    for (int f=0; f<n_flavors; ++f) {
      eps[f].resize(1);
      epssq[f].resize(1);
      eps[f][0] = parms["EPS_"+boost::lexical_cast<std::string>(f)];
      epssq[f][0] = parms["EPSSQ_"+boost::lexical_cast<std::string>(f)];
    }
    fourier_ptr.reset(new GFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., U,
                                              n_flavors, 1, densities, eps, epssq));
  }
  else if (parms.defined("CLUSTER_LOOP")){
    if(std::string(parms["CLUSTER_LOOP"])==std::string("DCA")){
      boost::shared_ptr<DCATransformer> clusterhandler;
      DCATransformer::generate_transformer(parms, clusterhandler);
      std::vector<std::vector<double> > eps_f(n_flavors, clusterhandler->epsav());
      std::vector<std::vector<double> > epssq_f(n_flavors, clusterhandler->epssqav());
      fourier_ptr.reset(new GFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., U,
                                                n_flavors, n_site, densities, eps_f, epssq_f));
    }else if(std::string(parms["CLUSTER_LOOP"])==std::string("CDMFT")){
      std::cout<<"using CDMFT cluster Fourier transformer"<<std::endl;
      std::cout<<"warning: G0 moments are not implemented for CDMFT"<<std::endl;
      std::vector<std::vector<double> > eps_f(n_flavors, std::vector<double>(n_site*n_site,0));
      std::vector<std::vector<double> > epssq_f(n_flavors, std::vector<double>(n_site*n_site,0));
      fourier_ptr.reset(new ClusterG0FourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2, U,
                                                        n_flavors, n_site, eps_f, epssq_f));
    }
    else throw std::logic_error("please figure out how to do fourier transforms in your cluster framwork if it is not DCA (3)");
  }
  else {
    std::cout << "using Bethe lattice fourier transform" << "\n";
    std::vector<std::vector<double> > eps(n_flavors);
    std::vector<std::vector<double> >epssq(n_flavors);
    double t = parms["t"];
    for (int f=0; f<n_flavors; ++f) {
      eps[f].resize(1);
      epssq[f].resize(1);
      eps[f][0] = 0.;
      epssq[f][0] = t*t;
    }
    fourier_ptr.reset(new GFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., U,
                                              n_flavors, 1, densities, eps, epssq));
  }
}



void FourierTransformer::generate_transformer_U(const alps::Parameters &parms,
                                                boost::shared_ptr<FourierTransformer> &fourier_ptr,
                                                const std::vector<double> &densities, const std::vector<double> &magnetization)
{
  int n_flavors = parms.value_or_default("FLAVORS", 2);
  int n_site = parms.value_or_default("SITES", 1);
  double U = parms["U"];
  if (parms.defined("GENERAL_FOURIER_TRANSFORMER")) {
    std::cout << "using general fourier transformer" << "\n";
    std::vector<double> eps(n_flavors);
    std::vector<double> epssq(n_flavors);
    for (int f=0; f<n_flavors; ++f) {
      eps[f] = parms["EPS_"+boost::lexical_cast<std::string>(f)];
      epssq[f] = parms["EPSSQ_"+boost::lexical_cast<std::string>(f)];
    }
    fourier_ptr.reset(new SimpleAntiFerroGFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., (double)parms["H"], U,
                                                             n_flavors, densities, magnetization,  eps, epssq));
  }
  else if (parms.defined("CLUSTER_LOOP")){
    if(std::string(parms["CLUSTER_LOOP"])==std::string("DCA")){
      std::cout << "using cluster fourier transformer" << "\n";
      boost::shared_ptr<DCATransformer> clusterhandler;
      DCATransformer::generate_transformer(parms, clusterhandler);
      std::vector<std::vector<double> > eps_f(n_flavors, clusterhandler->epsav());
      std::vector<std::vector<double> > epssq_f(n_flavors, clusterhandler->epssqav());
      fourier_ptr.reset(new ClusterAntiFerroGFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., (double)parms["H"], U,
                                                                n_flavors, n_site, densities, magnetization, eps_f, epssq_f));
    }else if(std::string(parms["CLUSTER_LOOP"])==std::string("CDMFT")){
      std::cout<<"using CDMFT cluster Fourier transformer"<<std::endl;
      std::cout<<"warning: G0 moments are not implemented for CDMFT"<<std::endl;
      std::vector<std::vector<double> > eps_f(n_flavors, std::vector<double>(n_site*n_site,0));
      std::vector<std::vector<double> > epssq_f(n_flavors, std::vector<double>(n_site*n_site,0));
      fourier_ptr.reset(new ClusterG0FourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2, U,
                                                        n_flavors, n_site, eps_f, epssq_f));
    }
    else throw std::logic_error("please figure out how to do fourier transforms in your cluster framwork if it is not DCA (1)");
  } else {
    std::cout << "using AFM Bethe lattice fourier transform" << "\n";
    std::vector<double> eps(n_flavors);
    std::vector<double> epssq(n_flavors);
    double t = parms["t"];
    for (int f=0; f<n_flavors; ++f) {
      eps[f] = 0.;
      epssq[f] = t*t;
    }
    fourier_ptr.reset(new SimpleAntiFerroGFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., (double)parms["H"], U,
                                                             n_flavors, densities, magnetization, eps, epssq));
  }
}




/*void FittingFourierTransformer::frequency_to_time_ft(itime_green_function_t &G_tau,
 const matsubara_green_function_t &G_omega) const
 {
 //fit c1_ to G_omega
 std::cout << "Fitting Fourier Transformer:" << std::endl;
 for (unsigned int f=0; f<G_omega.nflavor(); ++f) {
 for (unsigned int s=0; s<G_omega.nsite(); ++s) {
 std::complex<double> t1 = 0;
 std::complex<double> t2 = 0;
 for (unsigned int n=G_omega.nfreq()*3/4; n<G_omega.nfreq(); ++n) {
 std::complex<double> iw(0,(2*n+1)*M_PI/beta_);
 t1 += G_omega(n,s,s,f).imag()/iw;
 t2 += 1./(iw*iw);
 }
 std::complex<double> cc = t1/t2;
 c1_[f][s][s] = -cc.imag();
 std::cout << "c1(" << s << "," << s << "," << f << ")= " << c1_[f][s][s] << std::endl;
 }
 }
 FourierTransformer::frequency_to_time_ft(G_tau, G_omega);
 //calculating F(beta)
 unsigned int N_tau = G_tau.ntime()-1;
 unsigned int N_omega = G_omega.nfreq();
 unsigned int N_site = G_omega.nsite();
 matsubara_green_function_t G_omega_no_model(G_omega);
 double beta_minus_eps = beta_ - 1e-8;
 for(unsigned int f=0;f<G_omega.nflavor();++f){
 for (unsigned int s1=0; s1<N_site; ++s1){
 for (unsigned int s2=0; s2<N_site; ++s2) {
 for (unsigned int k=0; k<N_omega; k++) {
 std::complex<double> iw(0,(2*k+1)*M_PI/beta_);
 G_omega_no_model(k,s1,s2,f) -= f_omega(iw, c1_[f][s1][s2],c2_[f][s1][s2], c3_[f][s1][s2]);
 }
 G_tau(N_tau,s1,s2,f) = f_tau(beta_minus_eps, beta_, c1_[f][s1][s2], c2_[f][s1][s2], c3_[f][s1][s2]);
 for (unsigned int k=0; k<N_omega; k++) {
 double wt((2*k+1)*beta_minus_eps/beta_*M_PI);
 G_tau(N_tau,s1,s2,f) += 2/beta_*(cos(wt)*G_omega_no_model(k,s1,s2,f).real()+
 sin(wt)*G_omega_no_model(k,s1,s2,f).imag());
 }
 }
 }
 }
 }*/
