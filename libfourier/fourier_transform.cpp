/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2015 Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/

#include "alps/config.h" // needed to set up correct bindings

#include "boost/numeric/bindings/ublas.hpp"
#include "boost/numeric/bindings/lapack/driver/gesv.hpp"

#include "fourier_transform.h"

typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;

///same as the template function in fouriertransform.h, but not for the
//multiple vectors. Use green_function_t instead.
void FourierTransformer::frequency_to_time_ft(itime_green_function_t &G_tau,
                                              const matsubara_green_function_t &G_omega) const {
  if(G_tau.nflavor()!=G_omega.nflavor() || G_tau.nsite()!=G_omega.nsite()) throw std::logic_error("FT dimensions don't match.");
  unsigned int N_tau = G_tau.ntime()-1;
  unsigned int N_omega = G_omega.nfreq();
  unsigned int N_site = G_omega.nsite();
  matsubara_green_function_t G_omega_no_model(G_omega);

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
          if (!is_nonzero){
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
  std::cout<<"Foureir transform of nambu"<<std::endl;
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
      /*if((c1_[0][s1][s2]==0) && (c2_[0][s1][s2] == 0) && (c3_[0][s1][s2]==0)){  //nothing happening in this gf - this fct is completely zero.
        for (unsigned int i=0; i<=N_tau; i++) {
          G_tau(i,s1,s2,0,0)=0.;
          G_tau(i,s1,s2,0,1)=0.;
          G_tau(i,s1,s2,1,0)=0.;
          G_tau(i,s1,s2,1,1)=0.;
        }
      }*/
      //else
      {
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
  std::size_t Np1 = spline_matrix.size1();
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
  std::size_t Np1 = spline_matrix.size1();
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
  std::size_t Np1 = v.size();
  std::size_t N = Np1-1;
  std::size_t N_omega = v_omega.size();
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
  std::size_t Np1 = v.size();
  std::size_t N = Np1-1;
  std::size_t N_omega = v_omega.size();
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
  for (int flavor=0; flavor<G0_omega.nflavor(); ++flavor) {
    for (int k=0; k<G0_omega.nsite(); ++k) {
      std::cout << "append tail to self-energy with coefficients: "
      << " " << Sc0_[flavor][k][k]
      << " " << Sc1_[flavor][k][k]
      << " " << Sc2_[flavor][k][k] << std::endl;
      for (int freq=nfreq_measured; freq<G0_omega.nfreq(); ++freq) {
        std::complex<double> iw(0,(2*freq+1)*M_PI/beta_);
        std::complex<double> Sigma = Sc0_[flavor][k][k] + Sc1_[flavor][k][k]/iw + Sc2_[flavor][k][k]/(iw*iw);
        G_omega(freq, k, k, flavor) = 1./(1./G0_omega(freq, k, k, flavor) - Sigma);
      }
    }
  }
}


/*
void ClusterAntiFerroGFourierTransformer::append_tail(matsubara_green_function_t& G_omega,
                                                      const matsubara_green_function_t& G0_omega,
                                                      const unsigned int nfreq_measured) const
{
  matsubara_green_function_t selfenergy_tail(G0_omega.nfreq()-nfreq_measured, G0_omega.nsite(), G0_omega.nflavor(), G_omega.shape());
  selfenergy_tail.clear();
  matsubara_green_function_t G_omega_tail(selfenergy_tail);
  for (int flavor=0; flavor<G0_omega.nflavor(); ++flavor) {
    for (int freq=nfreq_measured; freq<G0_omega.nfreq(); ++freq) {
      std::complex<double> iw(0,(2*freq+1)*M_PI/beta_);
      for (int k=0; k<G0_omega.nsite(); k+=2) {
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
  for (int flavor=0; flavor<G0_omega.nflavor(); ++flavor) {
    for (int freq=0; freq<nfreq_measured; ++freq) {
      for (int k=0; k<G0_omega.nsite(); k+=2) {
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
  for (int flavor=0; flavor<G0_omega.nflavor(); ++flavor){
    for (int freq=nfreq_measured; freq<G0_omega.nfreq(); ++freq){
      for (int k=0; k<G0_omega.nsite(); k+=2){
        for (int p1=0; p1<=1; p1++){
          for (int p2=0; p2<=1; p2++){
            G_omega(freq, k+p1, k+p2, flavor) = G_omega_tail(freq-nfreq_measured, k+p1, k+p2, flavor);
          }
        }
      }
    }
  }
}*/
