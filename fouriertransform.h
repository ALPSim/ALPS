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
#ifndef FOURIER_TRANSFORM_H
#define FOURIER_TRANSFORM_H
#include <complex>
#include "types.h" // for the multiple_vector_type
#include "green_function.h" // for the multiple_vector_type
#include <alps/parameter.h>

#include <alps/alea.h>
#include <alps/alea/simpleobseval.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/vector.h>

typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;

void generate_spline_matrix(dense_matrix & spline_matrix);

void solve_second_derivatives(dense_matrix &spline_matrix, std::vector<double> &rhs);

//example:  V = vector_type
/*template <class V> 
void evaluate_second_derivatives(double dt, double c2, double c3, 
				 dense_matrix & spline_matrix, 
				 const V & g, V & second_derivatives) {
  
  // g, rhs and second_derivatives have dimension N+1
  unsigned int Np1 = spline_matrix.size1();
  double dtsq=dt*dt;
  
  // rhs is the vector containing the data of the curve y = g(tau), which allows to 
  // compute the vector of second derivatives y'' at times tau_n by evaluating
  // y'' = spline_matrix * rhs(y)
  //
  //                   (G''(0)+G''(beta))*(dt)^2                                
  //                   y0 - 2*y1 + y2
  //                   y1 - 2*y2 + y3
  // rhs = 6/(dt)^2 * ...
  //
  //                   yNp1-3 - 2*yNp1-2 + yNp1-1
  //                   (G'(0)+G'(beta))*(dt)^2 + y0 - y1 + yNp1-2 - yNp1-1     

  V rhs;
  rhs.resize(Np1);
  for(unsigned int i=0;i<Np1;++i){
    second_derivatives[i]=0;
  }
  rhs[0] = -c3*dtsq;  //think more about this!
  for (unsigned int i=1; i<Np1-1; i++) {
    rhs[i] = g[i-1]-2*g[i]+g[i+1];
  }
  rhs[Np1-1] = c2*dtsq + g[0] - g[1] + g[Np1-2] - g[Np1-1];
  solve_second_derivatives(spline_matrix, rhs);
  for(unsigned int i=0;i<Np1;++i){
    second_derivatives[i]=rhs[i]*6./dtsq;
  }
}*/


inline std::complex<double> f_omega(std::complex<double> iw, double c1, double c2, double c3) {
  std::complex<double> iwsq=iw*iw;
  return c1/iw + c2/(iwsq) + c3/(iw*iwsq);
}


inline double f_tau(double tau, double beta, double c1, double c2, double c3) {
  return -0.5*c1 + (c2*0.25)*(-beta+2.*tau) + (c3*0.25)*(beta*tau-tau*tau);
}


class FourierTransformer{
public:
  FourierTransformer(double beta, double mu, int n_flavor, int n_site, double epssq){
    beta_=beta;
    mu_=mu;
    c2_.resize(n_flavor);
    c3_.resize(n_flavor);
    c1_.resize(n_flavor);
    for(int f=0;f<n_flavor;++f){
      c2_[f].resize(n_site);
      c3_[f].resize(n_site);
      c1_[f].resize(n_site);
      for(int i=0;i<n_site;++i){
        c2_[f][i].resize(n_site);
        c3_[f][i].resize(n_site);
        c1_[f][i].resize(n_site);
        for(int j=0;j<n_site;++j){
          c1_[f][i][j]=i==j?1:0;
          //we need to rederive c2 and c3 for the most general case!
          c2_[f][i][j]=-mu;
          c3_[f][i][j]=epssq+mu*mu;
        }
      }
    }
  }
  FourierTransformer(double beta, double mu, double U, int n_flavor, int n_site, const std::vector<std::vector<std::vector<double> > > &densities, double epssq){
    beta_=beta;
    mu_=mu;
    c2_.resize(n_flavor);
    c3_.resize(n_flavor);
    c1_.resize(n_flavor);
    for(int f=0;f<n_flavor;++f){
      c2_[f].resize(n_site);
      c3_[f].resize(n_site);
      c1_[f].resize(n_site);
      for(int i=0;i<n_site;++i){
        c2_[f][i].resize(n_site);
        c3_[f][i].resize(n_site);
        c1_[f][i].resize(n_site);
        for(int j=0;j<n_site;++j){
          c1_[f][i][j]=i==j?1:0;
          //we need to rederive c2 and c3 for the most general case!
          c2_[f][i][j]=-mu+U*densities[f==1?0:1][i][j]; //formula 13
          double s = U*U*densities[f==1?0:1][i][j]*(1.-densities[f==1?0:1][i][j]);
          c3_[f][i][j]=epssq+c2_[f][i][j]*c2_[f][i][j] +s;
        }
      }
    }
  }
  FourierTransformer(const hifreq_moments &moments, const double &beta, const double &mu){
    c1_=moments.c1;
    c2_=moments.c2;
    c3_=moments.c3;
    beta_=beta;
    mu_=mu;
  }
  virtual ~FourierTransformer(){
  }
  virtual void forward_ft(const itime_green_function_t &G_tau, matsubara_green_function_t &G_omega) const;
  virtual void backward_ft(itime_green_function_t &G_tau, const matsubara_green_function_t &G_omega) const;
  virtual void backward_ft(green_function<alps::RealObsevaluator> &G_tau, 
  const green_function<alps::RealObsevaluator> &G_omega_real,
  const green_function<alps::RealObsevaluator> &G_omega_imag) const;
  static void generate_transformer(const alps::Parameters &parms,boost::shared_ptr<FourierTransformer> & fourier_ptr, double mu_shift=0);
  static void generate_transformer(const alps::Parameters &parms,boost::shared_ptr<FourierTransformer> & fourier_ptr, const hifreq_moments &moments);
  static void generate_transformer_U(const alps::Parameters &parms,boost::shared_ptr<FourierTransformer> & fourier_ptr,
  const std::vector<std::vector<std::vector<double> > > &densities);
  void setmu(double mu){ mu_=mu;}
protected:
  double beta_;
  double mu_;
  std::vector<std::vector<std::vector<double> > > c1_;
  std::vector<std::vector<std::vector<double> > > c2_;
  std::vector<std::vector<std::vector<double> > > c3_;
};


class BetheFourierTransformer:public FourierTransformer{
public:
  BetheFourierTransformer(double beta, double mu, double t, int n_flavor, int n_site)
    :FourierTransformer(beta, mu, n_flavor, n_site, t*t){}
  BetheFourierTransformer(double beta, double mu, double t, double U, const std::vector<std::vector<std::vector<double> > > &densities, int n_flavor, int n_site)
    :FourierTransformer(beta, mu, U, n_flavor, n_site, densities, t*t){
    std::cerr<<"using Bethe Fourier Transformer! "<<std::endl;
    std::cerr<<"mu=: "<<mu_<<std::endl;
    std::cerr<<"beta=: "<<beta_<<std::endl;
    std::cerr<<"U=: "<<U_<<std::endl;
  }
  virtual ~BetheFourierTransformer(){}
private:
  double t_;
  double U_;
};


class FFunctionFourierTransformer:public FourierTransformer{
public:
  FFunctionFourierTransformer(double beta, double mu, double epsilonsq_av, int n_flavor, int n_site)
    :FourierTransformer(beta,mu, n_flavor, n_site, 0){
    std::cout<<"FFourier Transformer: beta: "<<beta<<" mu: "<<mu<<"epsilonsq_av: "<<epsilonsq_av<<std::endl;
    epsilonsq_av_=epsilonsq_av; //this is the integral of the second moment of the dos: \int_-\infty^\infty e^2 rho(e) de. It is t^2 for semicircle...
    for(int f=0;f<n_flavor;++f){
      for(int i=0;i<n_site;++i){
        for(int j=0;j<n_site;++j){
          c1_[f][i][j]=epsilonsq_av;
          c2_[f][i][j]=0;//-(2*epsilonsq_av*mu+mu*mu*mu);
          c3_[f][i][j]=0;
        }
      }
    }
  }
  virtual void forward_ft(const itime_green_function_t &G_tau, matsubara_green_function_t &G_omega) const{
    std::cerr<<"implement 2nd derivatives for F function before using this!"<<std::endl;
    FourierTransformer::forward_ft(G_tau, G_omega);
  }
  virtual void backward_ft(itime_green_function_t &G_tau, const matsubara_green_function_t &G_omega) const{
    FourierTransformer::backward_ft(G_tau, G_omega);
    std::cout<<"correcting for F function"<<std::endl;
    for(unsigned int j=0;j<G_tau.nflavor();++j){
      G_tau(G_tau.ntime()-1,j)=-epsilonsq_av_-G_tau(0,j);
    }
  }
  virtual void backward_ft(green_function<alps::RealObsevaluator> &G_tau, 
			   const green_function<alps::RealObsevaluator> &G_omega_real,
			   const green_function<alps::RealObsevaluator> &G_omega_imag) const{

    std::cerr<<"backward ft: "<<&G_tau<<std::endl;
    std::cerr<<"backward ft: "<<&G_omega_real<<std::endl;
    std::cerr<<"backward ft: "<<&G_omega_imag<<std::endl;
    std::cerr<<"not implemented for obsevaluators. exiting."<<std::endl;
    abort();
  } //not implemented.

  virtual ~FFunctionFourierTransformer(){}
private:
  double epsilonsq_av_;
};



class FittingFourierTransformer : public FourierTransformer
{
  //a simple Fourier Transformer which fits c1_/w to the tail of the matsubara green function
public:
  
  FittingFourierTransformer(double beta, double mu, int n_flavor, int n_site)
    : FourierTransformer(beta, mu, n_flavor, n_site, 0)
  {
    std::cout<<"Fitting FFourier Transformer" << std::endl;
    for(int f=0;f<n_flavor;++f){
      for(int i=0;i<n_site;++i){
	for(int j=0;j<n_site;++j){
	  c1_[f][i][j]=0;
	  c2_[f][i][j]=0;
	  c3_[f][i][j]=0;
	}
      }
    }
  }
  
  virtual void forward_ft(const itime_green_function_t &G_tau, matsubara_green_function_t &G_omega) const
  {
    std::cerr<<"implement 2nd derivatives for F function before using this!"<<std::endl;
    FourierTransformer::forward_ft(G_tau, G_omega);
  }
  
  virtual void backward_ft(itime_green_function_t &G_tau, const matsubara_green_function_t &G_omega);
  
  virtual void backward_ft(green_function<alps::RealObsevaluator> &G_tau, 
			   const green_function<alps::RealObsevaluator> &G_omega_real,
			   const green_function<alps::RealObsevaluator> &G_omega_imag) const
  {
    
    std::cerr<<"backward ft: "<<&G_tau<<std::endl;
    std::cerr<<"backward ft: "<<&G_omega_real<<std::endl;
    std::cerr<<"backward ft: "<<&G_omega_imag<<std::endl;
    std::cerr<<"not implemented for obsevaluators. exiting."<<std::endl;
    abort();
  } //not implemented.
  
  virtual ~FittingFourierTransformer(){}
  
};




#endif
