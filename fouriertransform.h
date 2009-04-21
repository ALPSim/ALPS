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


class FourierTransformer
{
public:
  
  FourierTransformer(double beta, double mu, int n_flavor, int n_site)
  {
    init(beta, mu, n_flavor, n_site);
  }


  FourierTransformer(double beta, double mu, int n_flavor, int n_site, 
		     std::vector<std::vector<double> >& eps, std::vector<std::vector<double> >& epssq)
  {
    init(beta, mu, n_flavor, n_site);
    for(int f=0;f<n_flavor;++f){
      for(int i=0;i<n_site;++i){
	c1_[f][i][i] = 1.;
	c2_[f][i][i] = eps[f][i] - mu;
	c3_[f][i][i] = epssq[f][i] - 2*mu*eps[f][i] + mu*mu;
      }
    }
  }

  
  FourierTransformer(double beta, double mu, double U, int n_flavor, int n_site, 
		     const std::vector<double>& densities, 
		     std::vector<std::vector<double> >& eps, std::vector<std::vector<double> >& epssq)
  {
    init(beta, mu, n_flavor, n_site);
    for(int f=0;f<n_flavor;++f){
      int fbar = f%2==0 ? f+1 : f-1;
      std::cerr << "dens: " << densities[fbar] << std::endl;
      for(int i=0;i<n_site;++i){
	c1_[f][i][i] = 1.;
	c2_[f][i][i] = eps[f][i] - mu + U * densities[fbar]; //formula 13
	c3_[f][i][i] = epssq[f][i] - 2.*mu*eps[f][i] + mu*mu 
	  + 2.*U*densities[fbar]*(eps[f][i]-mu) + U*U*densities[fbar];
	Sc0_[f][i][i] = U * (densities[fbar]-0.5);
	Sc1_[f][i][i] = U*U * densities[fbar] * (1-densities[fbar]);
      }
    }
    /*std::ofstream forstr("for.dat");
    for(int f=0;f<n_flavor;++f){
      std::cout << "densities: " << densities[f] << std::endl; 
      for(int s1=0;s1<n_site;++s1){
	for(uint k=0; k<100; k++){
	  std::complex<double> iw(0,(2*k+1)*M_PI/beta_);
	  forstr << iw.imag() << "\t" << f << "\t" << s1 << "\t" 
		 << f_omega(iw, c1_[f][s1][s1],c2_[f][s1][s1], c3_[f][s1][s1]).real() << "\t"
		 << f_omega(iw, c1_[f][s1][s1],c2_[f][s1][s1], c3_[f][s1][s1]).imag()
		 << std::endl;
	}
      }
      }*/
  }

  virtual ~FourierTransformer() {}
  virtual void forward_ft(const itime_green_function_t &G_tau, matsubara_green_function_t &G_omega) const;
  virtual void backward_ft(itime_green_function_t &G_tau, const matsubara_green_function_t &G_omega) const;
  virtual void backward_ft(green_function<alps::RealObsevaluator> &G_tau, 
			   const green_function<alps::RealObsevaluator> &G_omega_real,
			   const green_function<alps::RealObsevaluator> &G_omega_imag) const;
  virtual void append_tail(matsubara_green_function_t& G_omega, const matsubara_green_function_t& G0_omega,
			   const int nfreq_measured) const;
  
  static void generate_transformer(const alps::Parameters &parms,
				   boost::shared_ptr<FourierTransformer> &fourier_ptr, 
				   std::vector<double> eps, std::vector<double> epssq);
  static void generate_transformer_U(const alps::Parameters &parms,
				     boost::shared_ptr<FourierTransformer> &fourier_ptr,
				     const std::vector<double> &densities,
				     std::vector<double> eps, std::vector<double> epssq);
  static void generate_transformer(const alps::Parameters &parms,
				   boost::shared_ptr<FourierTransformer> &fourier_ptr);
  static void generate_transformer_U(const alps::Parameters &parms,
				     boost::shared_ptr<FourierTransformer> &fourier_ptr,
				     const std::vector<double> &densities);
  //void setmu(double mu){ mu_=mu;}
  
protected:

  double beta_;
  double mu_;
  std::vector<std::vector<std::vector<double> > > c1_;
  std::vector<std::vector<std::vector<double> > > c2_;
  std::vector<std::vector<std::vector<double> > > c3_;
  std::vector<std::vector<std::vector<double> > > Sc0_; //coefficients for the self-energy
  std::vector<std::vector<std::vector<double> > > Sc1_;  

private:

  void init(double beta, double mu, int n_flavor, int n_site)
  {
    beta_=beta;
    mu_=mu;
    c1_.resize(n_flavor);
    c2_.resize(n_flavor);
    c3_.resize(n_flavor);
    Sc0_.resize(n_flavor);
    Sc1_.resize(n_flavor);
    for(int f=0;f<n_flavor;++f){
      c1_[f].resize(n_site);
      c2_[f].resize(n_site);
      c3_[f].resize(n_site);
      Sc0_[f].resize(n_site);
      Sc1_[f].resize(n_site);
      for(int i=0;i<n_site;++i){
        c1_[f][i].resize(n_site);
        c2_[f][i].resize(n_site);
        c3_[f][i].resize(n_site);
        Sc0_[f][i].resize(n_site);
        Sc1_[f][i].resize(n_site);
        for(int j=0;j<n_site;++j){
          c1_[f][i][j] = (i==j) ? 1. : 0.;
          c2_[f][i][j]= 0.;
          c3_[f][i][j]= 0.;
          Sc0_[f][i][j] = 0.;
          Sc1_[f][i][j]= 0.;
        }
      }
    } 
  }

};



class FFunctionFourierTransformer:public FourierTransformer
{
public:
  FFunctionFourierTransformer(double beta, double mu, double epsilonsq_av, int n_flavor, int n_site)
    :FourierTransformer(beta,mu, n_flavor, n_site){
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
    std::cout<<"implement 2nd derivatives for F function before using this!"<<std::endl;
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
    
    std::cout<<"backward ft: "<<&G_tau<<std::endl;
    std::cout<<"backward ft: "<<&G_omega_real<<std::endl;
    std::cout<<"backward ft: "<<&G_omega_imag<<std::endl;
    std::cout<<"not implemented for obsevaluators. exiting."<<std::endl;
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
    : FourierTransformer(beta, mu, n_flavor, n_site)
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
    std::cout<<"implement 2nd derivatives for F function before using this!"<<std::endl;
    FourierTransformer::forward_ft(G_tau, G_omega);
  }
  
  virtual void backward_ft(itime_green_function_t &G_tau, const matsubara_green_function_t &G_omega);
  
  virtual void backward_ft(green_function<alps::RealObsevaluator> &G_tau, 
			   const green_function<alps::RealObsevaluator> &G_omega_real,
			   const green_function<alps::RealObsevaluator> &G_omega_imag) const
  {
    
    std::cout<<"backward ft: "<<&G_tau<<std::endl;
    std::cout<<"backward ft: "<<&G_omega_real<<std::endl;
    std::cout<<"backward ft: "<<&G_omega_imag<<std::endl;
    std::cout<<"not implemented for obsevaluators. exiting."<<std::endl;
    abort();
  } //not implemented.
  
  virtual ~FittingFourierTransformer(){}
  
};




#endif
