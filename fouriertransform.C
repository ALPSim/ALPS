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
#include "cluster.h"
#include <valarray>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <alps/parameter.h>



/*void generate_spline_matrix(dense_matrix & spline_matrix) {
  
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
  //Emanuel avoid inversion of matrix
  spline_matrix=A;
}*/



void solve_second_derivatives(dense_matrix &spline_matrix, std::vector<double> &rhs){
  //avoid matrix inversion and instead solve the system of equations: S *v2 =
  //rhs for v2. The result will be found in rhs as we use the LAPACK call dgesv
  dense_matrix v(rhs.size(), 1);
  for(unsigned int i=0;i<rhs.size();++i){
    v(i,0)=rhs[i];
  }
  boost::numeric::bindings::lapack::gesv(spline_matrix, v);
  for(unsigned int i=0;i<rhs.size();++i){
    rhs[i]=v(i,0);
  }
}



void FourierTransformer::backward_ft(green_function<alps::RealObsevaluator> &G_tau, 
                                     const green_function<alps::RealObsevaluator> &G_omega_real,
                                     const green_function<alps::RealObsevaluator> &G_omega_imag) const
{
  assert(G_tau.nflavor()==G_omega_real.nflavor() && G_tau.nsite()==G_omega_real.nsite());
  assert(G_tau.nflavor()==G_omega_imag.nflavor() && G_tau.nsite()==G_omega_imag.nsite());
  unsigned int N_tau = G_tau.ntime()-1;
  unsigned int N_omega = G_omega_real.nfreq();
  unsigned int N_site = G_omega_real.nsite();
  unsigned int N_flavor = G_omega_real.nflavor();
  green_function<alps::RealObsevaluator> G_omega_no_model_real(G_omega_real);
  green_function<alps::RealObsevaluator> G_omega_no_model_imag(G_omega_imag);
  double dt = beta_/N_tau;
  for(unsigned int f=0;f<N_flavor;++f){
    for (unsigned int s1=0; s1<N_site; ++s1){
      for (unsigned int s2=0; s2<N_site; ++s2) {
        if (s1==s2) {
          for (unsigned int k=0; k<N_omega; k++) {
            std::complex<double> iw(0,(2*k+1)*M_PI/beta_);
            G_omega_no_model_real(k,s1,s2,f) -= f_omega(iw, c1_[f][s1][s2],c2_[f][s1][s2], c3_[f][s1][s2]).real();
            G_omega_no_model_imag(k,s1,s2,f) -= f_omega(iw, c1_[f][s1][s2],c2_[f][s1][s2], c3_[f][s1][s2]).imag();
          } 
        }
        for (unsigned int i=0; i<N_tau; i++) {
          double wt = i*M_PI/N_tau;
          G_tau(i,s1,s2,f) = 2/beta_*(cos(wt)*G_omega_no_model_real(0,s1,s2,f)+
                                      sin(wt)*G_omega_no_model_imag(0,s1,s2,f));
          for (unsigned int k=1; k<N_omega; k++) {
            wt = (2*k+1)*i*M_PI/N_tau;
            G_tau(i,s1,s2,f) += 2/beta_*(cos(wt)*G_omega_no_model_real(k,s1,s2,f)+
                                         sin(wt)*G_omega_no_model_imag(k,s1,s2,f));
          }
          G_tau(i,s1,s2,f) += f_tau(i*dt, beta_, c1_[f][s1][s2], c2_[f][s1][s2], c3_[f][s1][s2]); 
        }
        G_tau(N_tau,s1,s2,f)= -G_tau(0,s1,s2,f);
        G_tau(N_tau,s1,s2,f)-= s1==s2 ? -1. : 0.;
      }
    }
  }
}



///same as the template function in fouriertransform.h, but not for the
//multiple vectors. Use green_function_t instead.
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
  //std::cout<<"Fourier: G omega nomodel: "<<G_omega_no_model<<std::endl;
  //std::cout<<"this was generated for c1: "<<c1_[0][0][0]<<"c2: "<<c2_[0][0][0]<<std::endl;

  /*std::ofstream backstr("back.dat");
  G_omega.write("gtest");
  test.write("ftest");
  int sdel = 1;//N_site > 1 ? 2 : 1;
  int pmax = 0;//N_site > 1 ? 1 : 0;
  for(unsigned int f=0;f<G_omega.nflavor();++f){
    for (unsigned int s=0; s<N_site; s+=sdel){
      for (int p1=0; p1<=pmax; ++p1) {
        for (int p2=0; p2<=pmax; ++p2) {
          for(uint k=0; k<N_omega; k++){
            std::complex<double> iw(0,(2*k+1)*M_PI/beta_);
            backstr << iw.imag() << "\t" << f << "\t" << s+p1 << "\t" << s+p2 << "\t"
            << G_omega(k,s+p1,s+p2,f).real() << "\t" << G_omega(k,s+p1,s+p2,f).imag() << "\t" 
            << f_omega(iw, c1_[f][s+p1][s+p2],c2_[f][s+p1][s+p2], c3_[f][s+p1][s+p2]).real() << "\t"
            << f_omega(iw, c1_[f][s+p1][s+p2],c2_[f][s+p1][s+p2], c3_[f][s+p1][s+p2]).imag()
            << std::endl;
          }
        }
      }
    }
  }*/
}
  
  



/*void FourierTransformer::forward_ft(const itime_green_function_t & G_tau, matsubara_green_function_t & G_omega) const 
{
  assert(G_tau.nflavor()==G_omega.nflavor() && G_tau.nsite()==G_omega.nsite());
  uint N=G_tau.ntime()-1;
  uint N_site=G_tau.nsite();
  uint n_oversampling_t=10; //not sure that helps... check again!
  uint n_oversampling_w=1; //not sure that helps... check again!
  uint N_large_t=(G_tau.ntime()-1)*n_oversampling_t+1;
  uint N_large_w=(G_tau.ntime()-1)*n_oversampling_w+1;
  double dt=beta_/(N_large_t-1);
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
	  double phi=w_n*i*dt;
	  G_omega_vector[n] += std::complex<double>(cos(phi), sin(phi))*
	  ((G_tau_vector[i]-f_tau(dt*i,beta_, c1_[flavor][s1][s2], c2_[flavor][s1][s2], c3_[flavor][s1][s2]))*dt);
	  }
	  }
        std::cout << "running forward ft with coefficients: "
	<< c1_[flavor][s1][s2] << " " << c2_[flavor][s1][s2] << " "
	<< c3_[flavor][s1][s2] << std::endl;
	for(uint n=0;n<G_omega.nfreq();++n){
	G_omega(n,s1,s2,flavor)=n<N_large_w?G_omega_vector[n]:0.;
	std::complex<double> iw(0,(2*n+1)*M_PI/beta_);
	  G_omega(n,s1,s2,flavor) += f_omega(iw, c1_[flavor][s1][s2], c2_[flavor][s1][s2], c3_[flavor][s1][s2]);
	}
      }
    }
  }
}


 
void FourierTransformer::forward_ft(const itime_green_function_t & G_tau, matsubara_green_function_t & G_omega) const {
  assert(G_tau.nflavor()==G_omega.nflavor() && G_tau.nsite()==G_omega.nsite());
  unsigned int N=G_tau.ntime()-1;
  unsigned int n_omega=G_omega.nfreq();
  unsigned int N_site=G_tau.nsite();
  double dt=beta_/N;
  //we implicitely assume that we have as many imag time slices as matsubara frequencies. 
  //If this is not the case: interpolate,
  double G_tau_vector[N];
  std::valarray<std::complex<double> > G_omega_vector(n_omega);
  std::valarray<std::complex<double> > G_omega_vector_nosub(n_omega);
  //assert(c1_==1); //otherwise change evaluate_second_derivatives!
  for(unsigned int flavor=0;flavor<G_tau.nflavor();++flavor){
    for (unsigned int s1=0; s1<N_site; ++s1){
      for (unsigned int s2=0; s2<N_site; ++s2) {
        for(unsigned int i=0;i<N;++i){
          G_tau_vector[i]=G_tau(i,s1,s2,flavor);
        }
        for (unsigned int n=0; n<n_omega; n++) {
          double w_n=(2*n+1)*M_PI/beta_;
          G_omega_vector[n]=0;
          G_omega_vector_nosub[n]=0;
         for(unsigned int i=0;i<N;++i){
            double phi=w_n*i*dt;
	    G_omega_vector[n] += std::complex<double>(cos(phi), sin(phi))*
	      ((G_tau_vector[i]-f_tau(dt*i,beta_, c1_[flavor][s1][s2], c2_[flavor][s1][s2], c3_[flavor][s1][s2]))*dt);
	    G_omega_vector_nosub[n] += std::complex<double>(cos(phi), sin(phi))*G_tau_vector[i];
            //  G_omega_vector[n] += dt*std::exp((std::complex<double>(0,1))*(w_n*i*dt))*G_tau_vector[i];
          }
          if(n>n_omega/2)
            G_omega_vector[n]=0; //no more information for these frequencies.
        }
        std::cout<<"doing forward ft for nfreq= "<<G_omega.nfreq()<<std::endl;
        std::cout<<"c1: "<<c1_[flavor][s1][s2]<<" c2: "<<c2_[flavor][s1][s2]<<" c3: "<<c3_[flavor][s1][s2]<<std::endl;
        for(unsigned int n=0;n<G_omega.nfreq();++n){
          G_omega(n,s1,s2,flavor)=G_omega_vector[n];
	  std::complex<double> iw(0,(2*n+1)*M_PI/beta_);
	  G_omega(n,s1,s2,flavor) += f_omega(iw, c1_[flavor][s1][s2],c2_[flavor][s1][s2], c3_[flavor][s1][s2]);
        }
      }
    }
  }
  //for(unsigned int i=0;i<N;++i){
    //std::cerr << i*dt << "\t" 
	//      << G_tau_vector[i] << "\t"
	  //    << f_tau(dt*i,beta_, c1_[0][0][0], c2_[0][0][0], c3_[0][0][0]) << std::endl;
	    //  }
}*/



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
  boost::numeric::bindings::lapack::gesv(A, spline_matrix);
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
  /*for(int i=0;i<Np1;++i){
    g[i]-=f_tau(dt*i, BETA, c1g, c2g, c3g);
  }*/
  //rhs[0] = 10000*pow(1./100,2)*10000; // const2*dtau^2*const2 for G^-1=-const2*G; this is the condition if we do not remove the parabolic part !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //std::cout<<g[0]<<" "<<g[1]<<" "<<g[Np1-1]<<" "<<g[Np1-2]<<std::endl;
  rhs[0] = -c3g; //G''(0)+G''(beta)=-c3
  for (int i=1; i<Np1-1; i++) {
    rhs[i] = (g[i-1]-2*g[i]+g[i+1])/dt;
  }
  //rhs[Np1-1] = g[0] - g[1] + g[Np1-2] - g[Np1-1];
  rhs[Np1-1] = c2g -1./dt*(-g[0] + g[1] -g[Np1-2] + g[Np1-1]);

  
  //rhs[Np1-1] = 0; // natural splines
              
  for (int i=0; i<Np1; i++) {
    second_derivatives[i]=0;
    for (int j=0; j<Np1; j++) {
      second_derivatives[i] += spline_matrix(i,j)*rhs[j];
    }
  //  std::cout << "2der " << i << " " << second_derivatives[i] << "  rhs " << rhs[i] << " g "<<g[i]<<" f "<<f_tau(dt*i,BETA,c1g,c2g,c3g)<<"\n";
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
  
  for(uint f=0;f<gtau.nflavor();++f){
    for(uint p=0;p<gtau.nsite();++p){
      for(uint q=0;q<gtau.nsite();++q){
        for(int tau=0;tau<Np1;++tau){
          v[tau]=gtau(tau,p,q,f);
        }
        /*std::cout << "running forward ft with coefficients: "
		  << " " << c1_[f][p][q]
		  << " " << c2_[f][p][q]
		  << " " << c3_[f][p][q] << std::endl;*/
        //std::cout<<"c1 is: "<<c1_[f][p][q]<<" computed: "<<-v[0]-v[N]<<std::endl;
        //std::cout<<"c2 is: "<<c2_[f][p][q]<<" computed: "<<(v[1]-v[0]+v[N]-v[N-1])/dt<<std::endl;
        //std::cout<<"c3 is: "<<c3_[f][p][q]<<" computed: "<<-(v[2]-2*v[1]+v[0]+v[N-2]-2*v[N-1]+v[N])/(dt*dt)<<std::endl;

        //  *(const_cast<double*> (&(c1_[f][p][q])))=1;
        //  *(const_cast<double*> (&(c2_[f][p][q])))*=-1; //convention is all wrong...
        //  *(const_cast<double*> (&(c3_[f][p][q])))=0;

        //std::cout<<"c1 new is: "<<c1_[f][p][q]<<" computed: "<<-v[0]-v[N]<<std::endl;
        //std::cout<<"c2 new is: "<<c2_[f][p][q]<<" computed: "<<(v[1]-v[0]+v[N]-v[N-1])/dt<<std::endl;
        //std::cout<<"c3 new is: "<<c3_[f][p][q]<<" computed: "<<-(v[2]-2*v[1]+v[0]+v[N-2]-2*v[N-1]+v[N])/(dt*dt)<<std::endl;
       
        /*bool use_computed=true;
        if(use_computed){
          std::cout<<"using computed, not analytic values for c1g, c2g, c3g"<<std::endl;
          *(const_cast<double*> (&(c1_[f][p][q])))=-v[0]-v[N];
          *(const_cast<double*> (&(c2_[f][p][q])))=(v[1]-v[0]+v[N]-v[N-1])/dt;
          *(const_cast<double*> (&(c3_[f][p][q])))=-(v[2]-2*v[1]+v[0]+v[N-2]-2*v[N-1]+v[N])/(dt*dt);
	  }*/
        
        dense_matrix spline_matrix(Np1, Np1);
        generate_spline_matrix(spline_matrix, dt);
        // matrix containing the second derivatives y'' of interpolated y=v[tau] at points tau_n 
        std::vector<double> v2(Np1, 0); 
        evaluate_second_derivatives(dt/*,beta_*/, spline_matrix, v, v2, c1_[f][p][q], c2_[f][p][q], c3_[f][p][q]);
        v_omega.assign(N_omega, 0);
        //std::vector<std::complex<double> > v_omega_naive(N_omega);
        //std::vector<std::complex<double> > v_omega_third(N_omega);
        //v_omega_naive.assign(N_omega, 0);
        //v_omega_third.assign(N_omega, 0);
        
        //std::cout << "beta " << beta_ << "\n";
        
        for (int k=0; k<N_omega; k++) {
          std::complex<double> iw(0, M_PI*(2*k+1)/beta_);
          for (int n=1; n<N; n++) {
            v_omega[k] += exp(iw*(n*dt))*(v2[n+1]-2*v2[n]+v2[n-1]); //partial integration, four times. Then approximate the fourth derivative by finite differences
            //v_omega_naive[k]+=std::exp(iw*(n*dt))*v[n];
                                                                    // v_omega_third[k]+=std::exp(iw*(n*dt))*(-v2[n]+v2[n+1])/(iw*iw*iw);
          }
          v_omega[k] += (v2[1] - v2[0] + v2[N] - v2[N-1]); //that's the third derivative, on the boundary
          v_omega[k] *= 1./(dt*iw*iw*iw*iw);
          //v_omega_naive[k]+=v[0]-v[N];
          //v_omega_third[k]+=(-v2[0]+v2[1]); //first slice
          //v_omega_third[k]=-v_omega_third[k]+f_omega(iw, c1g[f], c2g[f], c3g[f]);
          //std::cout<<k<<" "<<v_omega_naive[k]<<" "<<v_omega[k]/*<<" "<<v_omega_third[k]*/<<" "<<f_omega(iw, c1_[f][p][q],c2_[f][p][q],c3_[f][p][q])<<std::endl;
          v_omega[k] += f_omega(iw, c1_[f][p][q], c2_[f][p][q], c3_[f][p][q]); //that's the boundary terms of the first, second, and third partial integration.
                                                                               //std::cout<<"derivative at boundary: "<<-v2[1] + v2[0] + v2[N] - v2[N-1]<<" divisor: "<< 1./(dt*iw*iw*iw*iw)<<std::endl;
          //gomega(k,f)=std::conj(v_omega[k]); //that's the proper convention for the self consistency loop in my DCA program.
          gomega(k,p,q,f)=v_omega[k]; //that's the proper convention for the self consistency loop here.
        }
        
        /*for (int k=0; k<N_omega; k++) {
          std::complex<double> fw=f_omega(std::complex<double>(0, (2.*k+1)*M_PI/beta_),c1_[f][p][q],c2_[f][p][q],c3_[f][p][q]);
	  std::cout << "v_omega " << k << " " << M_PI*(2*k+1) << " " << std::real(v_omega[k]) << "  " << std::imag(v_omega[k]) << " "<<fw.real()<<" "<<fw.imag()<<" "<<(v_omega[k]-fw).real()<<" "<<(v_omega[k]-fw).imag()<<" naive: "<<v_omega_naive[k].real()<<" "<<v_omega_naive[k].imag()<<std::endl;
        }*/
        
        /*static std::ofstream g_interpol("/tmp/g_interpol.dat");                      
        for (int j=0; j<Np1-1; j++) {
          // define parameters of spline fit y(tau) = a + b*(tau-tau_j) + c*(tau-tau_j)^2 + d*(tau-tau_j)^3
          double a = v[j];
          double b = (v[j+1]-v[j])/dt - dt/6.*(v2[j+1]+2.*v2[j]);
          double c = v2[j]/2.;
          double d = (v2[j+1]-v2[j])/(6.*dt);
          // plot interpolated function
          for (int n=0; n<10; n++) {
            double x = n/10.*dt;
            g_interpol<<j+n/10.<< " "<<a + b*x + c*pow(x,2) + d*pow(x,3)<<" "<<f_tau(x, beta_, c1_[f][p][q], c2_[f][p][q], c3_[f][p][q])<<std::endl;
          }
	  }
	  g_interpol<<std::endl;
	  static std::ofstream second_der("/tmp/2nd_der.dat");
	  for (int j=0; j<Np1; j++) {
          second_der<<j*dt<<" "<<v[j]<<" "<<v2[j]<<std::endl;
        }
        second_der<<std::endl;*/
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
  /*std::ofstream self("self.dat");
  for (spin_t flavor=0; flavor<G0_omega.nflavor(); ++flavor) {
    for (site_t k=0; k<G0_omega.nsite(); ++k) {
      for (frequency_t freq=0; freq<nfreq_measured; ++freq) {
	std::complex<double> Sigma = 1./G0_omega(freq, k, k, flavor) - 1./G_omega(freq, k, k, flavor);
	self << (2*freq+1)*M_PI/beta_ << " " << flavor << " " << k << " " 
	     << Sigma.real() << " " << Sigma.imag() << "\n";
      }
      for (frequency_t freq=nfreq_measured; freq<G0_omega.nfreq(); ++freq) {
	std::complex<double> iw(0,(2*freq+1)*M_PI/beta_);
	std::complex<double> Sigma = Sc0_[flavor][k][k] + Sc1_[flavor][k][k]/iw + Sc2_[flavor][k][k]/(iw*iw);
	self << (2*freq+1)*M_PI/beta_ << " " << flavor << " " << k << " " 
	     << Sigma.real() << " " << Sigma.imag() << "\n";
      }
    }
  }*/
}



void ClusterAntiFerroGFourierTransformer::append_tail(matsubara_green_function_t& G_omega, 
						      const matsubara_green_function_t& G0_omega,
						      const int nfreq_measured) const
{
  matsubara_green_function_t selfenergy_tail(G0_omega.nfreq()-nfreq_measured, G0_omega.nsite(), G0_omega.nflavor());
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
  matsubara_green_function_t selfenergy_measured(nfreq_measured, G0_omega.nsite(), G0_omega.nflavor());
  matsubara_green_function_t G_omega_measured(nfreq_measured, G0_omega.nsite(), G0_omega.nflavor());
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
  /*std::ofstream mappend("mappend.dat");
  for (site_t k=0; k<G0_omega.nsite(); k+=2) {
    for (int p1=0; p1<=1; p1++) {
      for (int p2=0; p2<=1; p2++) {
	for (spin_t flavor=0; flavor<G0_omega.nflavor(); ++flavor) {
	  for (frequency_t freq=0; freq<nfreq_measured; ++freq) {
	    std::complex<double> iw(0,(2*freq+1)*M_PI/beta_);
	    mappend << iw.imag() << "\t" << flavor << "\t" << k+p1 << "\t" << k+p2 << "\t"
		    << selfenergy_measured(freq, k+p1, k+p2, flavor).real() << "\t"
		    << selfenergy_measured(freq, k+p1, k+p2, flavor).imag() << std::endl;
	  }
	}
      }
    }
    }
  std::ofstream tappend("tappend.dat");
  for (site_t k=0; k<G0_omega.nsite(); k+=2) {
    for (int p1=0; p1<=1; p1++) {
      for (int p2=0; p2<=1; p2++) {
	for (spin_t flavor=0; flavor<G0_omega.nflavor(); ++flavor) {
	  for (frequency_t freq=nfreq_measured; freq<G0_omega.nfreq(); ++freq) {
	    std::complex<double> iw(0,(2*freq+1)*M_PI/beta_);
	    tappend << iw.imag() << "\t" << flavor << "\t" << k+p1 << "\t" << k+p2 << "\t"
		    << selfenergy_tail(freq-nfreq_measured, k+p1, k+p2, flavor).real() << "\t"
		    << selfenergy_tail(freq-nfreq_measured, k+p1, k+p2, flavor).imag() << std::endl;
	  }
	}
      }
    }
    }*/
  G_omega_tail = AntiFerroDCATransformer::green_function_inverse(G_omega_tail);
  G_omega_tail = AntiFerroDCATransformer::green_function_subtract(G_omega_tail, selfenergy_tail);
  G_omega_tail = AntiFerroDCATransformer::green_function_inverse(G_omega_tail);
  for (spin_t flavor=0; flavor<G0_omega.nflavor(); ++flavor) 
    for (frequency_t freq=nfreq_measured; freq<G0_omega.nfreq(); ++freq) 
      for (site_t k=0; k<G0_omega.nsite(); k+=2) 
	for (int p1=0; p1<=1; p1++) 
	  for (int p2=0; p2<=1; p2++) 
	    G_omega(freq, k+p1, k+p2, flavor) = G_omega_tail(freq-nfreq_measured, k+p1, k+p2, flavor);
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
    if(!(std::string(parms["CLUSTER_LOOP"])==std::string("DCA"))) throw std::logic_error("please figure out how to do fourier transforms in your cluster framwork if it is not DCA");
    std::cout << "using cluster fourier transformer" << "\n";
    boost::shared_ptr<DCATransformer> clusterhandler;
    DCATransformer::generate_transformer(parms, clusterhandler);
    std::vector<std::vector<double> > eps_f(n_flavors, clusterhandler->epsav());
    std::vector<std::vector<double> > epssq_f(n_flavors, clusterhandler->epssqav());
    fourier_ptr.reset(new ClusterG0FourierTransformer((double)parms["BETA"], (double)parms["MU"], (double)parms["H"], 
						      n_flavors, n_site, eps_f, epssq_f));    
  } 
  else {
    std::cout << "using Bethe lattice fourier transformer" << "\n";
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
    if(!(std::string(parms["CLUSTER_LOOP"])==std::string("DCA"))) throw std::logic_error("please figure out how to do fourier transforms in your cluster framwork if it is not DCA");
    std::cout << "using cluster fourier transformer" << "\n";
    boost::shared_ptr<DCATransformer> clusterhandler;
    DCATransformer::generate_transformer(parms, clusterhandler);
    std::vector<std::vector<double> > eps_f(n_flavors, clusterhandler->epsav());
    std::vector<std::vector<double> > epssq_f(n_flavors, clusterhandler->epssqav());
    fourier_ptr.reset(new GFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., U, 
					      n_flavors, n_site, densities, eps_f, epssq_f));    
  }
  else {
    std::cout << "using Bethe lattice fourier transformer" << "\n";
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
    if(!(std::string(parms["CLUSTER_LOOP"])==std::string("DCA"))) throw std::logic_error("please figure out how to do fourier transforms in your cluster framwork if it is not DCA");
    std::cout << "using cluster fourier transformer" << "\n";
    boost::shared_ptr<DCATransformer> clusterhandler;
    DCATransformer::generate_transformer(parms, clusterhandler);
    std::vector<std::vector<double> > eps_f(n_flavors, clusterhandler->epsav());
    std::vector<std::vector<double> > epssq_f(n_flavors, clusterhandler->epssqav());
    fourier_ptr.reset(new ClusterAntiFerroGFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., (double)parms["H"], U, 
							      n_flavors, n_site, densities, magnetization, eps_f, epssq_f));
  } else {
    std::cout << "using Bethe lattice fourier transformer" << "\n";
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




void FittingFourierTransformer::backward_ft(itime_green_function_t &G_tau, 
					    const matsubara_green_function_t &G_omega) 
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
  FourierTransformer::backward_ft(G_tau, G_omega);
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
}
