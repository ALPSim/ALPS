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
//#include "hirschfyesim.h"
#include <complex>
#include "types.h" // for the multiple_vector_type
#include "green_function.h" // for the multiple_vector_type
#include <alps/parameter.h>
typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;
void generate_spline_matrix(dense_matrix & spline_matrix);
void solve_second_derivatives(dense_matrix &spline_matrix, std::vector<double> &rhs);
//example:  V = vector_type
template <class V> void evaluate_second_derivatives(double dt, double c2, double c3, dense_matrix & spline_matrix, const V & g, V & second_derivatives) {

  // g, rhs and second_derivatives have dimension N+1
  uint Np1 = spline_matrix.size1();
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
  for(uint i=0;i<Np1;++i){
    second_derivatives[i]=0;
  }
  rhs[0] = -c3*dtsq;  //think more about this!
  for (uint i=1; i<Np1-1; i++) {
    rhs[i] = g[i-1]-2*g[i]+g[i+1];
  }
  rhs[Np1-1] = c2*dtsq + g[0] - g[1] + g[Np1-2] - g[Np1-1];
  //std::cout<<"rhs: "<<rhs<<std::endl;
  solve_second_derivatives(spline_matrix, rhs);
  //std::cout<<"rhs: "<<rhs<<std::endl;
  for(uint i=0;i<Np1;++i){
    second_derivatives[i]=rhs[i]*6./dtsq;
  }
  /*
  for (int i=0; i<Np1; i++) {
    for (int j=0; j<Np1; j++) {
      second_derivatives[i] += spline_matrix(i,j)*rhs[j];
	}
	second_derivatives[i] *= 6./dtsq;
  }*/
}

inline std::complex<double> f_omega(std::complex<double> iw, double c1, double c2, double c3) {
  std::complex<double> iwsq=iw*iw;
	return c1/iw + c2/(iwsq) + c3/(iw*iwsq);
}

inline double f_tau(double tau, double beta, double c1, double c2, double c3) {
  return -0.5*c1 + (c2*0.25)*(-beta+2.*tau) + (c3*0.25)*(beta*tau-tau*tau);
}
class FourierTransformer{
  public:
    FourierTransformer(double beta, double mu){
      beta_=beta;
      mu_=mu;
      c1_=1.;
      c2_=mu/*/beta*/;
      c3_=0.;
    }
    virtual ~FourierTransformer(){}
    virtual void forward_ft(const itime_green_function_t &G_tau, matsubara_green_function_t &G_omega) const;
    virtual void backward_ft(itime_green_function_t &G_tau, const matsubara_green_function_t &G_omega) const;
    static void generate_transformer(const alps::Parameters &parms,boost::shared_ptr<FourierTransformer> & fourier_ptr);
    static void generate_transformer_U(const alps::Parameters &parms,boost::shared_ptr<FourierTransformer> & fourier_ptr);
  protected:
    double beta_;
    double mu_;
    double c1_;
    double c2_;
    double c3_;
};

class BetheFourierTransformer:public FourierTransformer{
  public:
  BetheFourierTransformer(double beta, double mu, double t):FourierTransformer(beta, mu){
    t_=t;
    c3_=(mu_*mu_+t_*t_)/*/(beta*beta)*/;
    U_=0;
  }
  BetheFourierTransformer(double beta, double mu, double t, double U):FourierTransformer(beta, mu){
    t_=t;
    U_=U;
    c3_=(mu_*mu_+t_*t_+ U_*U_/4.)/*/(beta*beta)*/; //the U-trick: for the GF (not bgf)
    std::cerr<<"using Bethe Fourier Transformer! "<<std::endl;
    std::cerr<<"mu=: "<<mu_<<std::endl;
    std::cerr<<"beta=: "<<beta_<<std::endl;
    std::cerr<<"U=: "<<U_<<std::endl;
    std::cerr<<"c3=: "<<c3_<<std::endl;
  }
  virtual ~BetheFourierTransformer(){}
  private:
  double t_;
  double U_;
};
class FFunctionFourierTransformer:public FourierTransformer{
  public:
    FFunctionFourierTransformer(double beta, double mu, double t):FourierTransformer(beta,mu){
      std::cout<<"FFourier Transformer: beta: "<<beta<<" mu: "<<mu<<"t: "<<t<<std::endl;
      t_=t;
      c1_=t*t;
      c2_=-(2*t*t*mu+mu*mu*mu);
      c3_=0;//t*t*t*t-t*t*mu*mu-mu*mu*mu*mu; //recheck this term!
    }
    virtual void forward_ft(const itime_green_function_t &G_tau, matsubara_green_function_t &G_omega) const{
      std::cerr<<"implement 2nd derivatives for F function before using this!"<<std::endl;
      FourierTransformer::forward_ft(G_tau, G_omega);
    }
    virtual void backward_ft(itime_green_function_t &G_tau, const matsubara_green_function_t &G_omega) const{
      FourierTransformer::backward_ft(G_tau, G_omega);
      std::cout<<"correcting for F function"<<std::endl;
      for(uint j=0;j<G_tau.nflavor();++j){
        G_tau(G_tau.ntime()-1,j)=-t_*t_-G_tau(0,j);
        /*G_tau(0,j)=-t_*t_*(0.5+0.25*mu_);
        G_tau(G_tau.ntime()-1,j)=-t_*t_-G_tau(0,j);*/
      }
    }
    virtual ~FFunctionFourierTransformer(){}
  private:
    double t_;
};
#endif
