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
#include "boost/shared_ptr.hpp"
#include "types.h" 
#include "green_function.h" 
#include "nambu_green_function.h" 
#include "alps/ngs/params.hpp"


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
  
  FourierTransformer(const double beta, const int n_flavor, const int n_site)
  {
    beta_=beta;
    //mu_=mu;
    c1_.resize(n_flavor);
    c2_.resize(n_flavor);
    c3_.resize(n_flavor);
    Sc0_.resize(n_flavor);
    Sc1_.resize(n_flavor);
    Sc2_.resize(n_flavor);
    for(int f=0;f<n_flavor;++f){
      c1_[f].resize(n_site);
      c2_[f].resize(n_site);
      c3_[f].resize(n_site);
      Sc0_[f].resize(n_site);
      Sc1_[f].resize(n_site);
      Sc2_[f].resize(n_site);
      for(int i=0;i<n_site;++i){
        c1_[f][i].resize(n_site);
        c2_[f][i].resize(n_site);
        c3_[f][i].resize(n_site);
        Sc0_[f][i].resize(n_site);
        Sc1_[f][i].resize(n_site);
        Sc2_[f][i].resize(n_site);
        for(int j=0;j<n_site;++j){
          c1_[f][i][j] = (i==j) ? 1. : 0.;
          c2_[f][i][j]= 0.;
          c3_[f][i][j]= 0.;
          Sc0_[f][i][j] = 0.;
          Sc1_[f][i][j]= 0.;
          Sc2_[f][i][j]= 0.;
        }
      }
    } 
  }
  
  
  virtual ~FourierTransformer() {}
  //standard forward and backward transformer
  virtual void time_to_frequency_ft(const itime_green_function_t &G_tau, matsubara_green_function_t &G_omega) const;
  virtual void frequency_to_time_ft(itime_green_function_t &G_tau, const matsubara_green_function_t &G_omega) const;
  //NAMBU forward and backward transformer
  virtual void time_to_frequency_ft(const nambu_itime_green_function_t &G_tau, nambu_matsubara_green_function_t &G_omega) const;
  virtual void frequency_to_time_ft(nambu_itime_green_function_t &G_tau, const nambu_matsubara_green_function_t &G_omega) const;
  /*  virtual void frequency_to_time_ft(green_function<alps::RealObsevaluator> &G_tau, 
   const green_function<alps::RealObsevaluator> &G_omega_real,
   const green_function<alps::RealObsevaluator> &G_omega_imag) const;*/
  virtual void append_tail(matsubara_green_function_t& G_omega, const matsubara_green_function_t& G0_omega,
                           const unsigned int nfreq_measured) const;
  
  static void generate_transformer(const alps::params &parms,
                                   boost::shared_ptr<FourierTransformer> &fourier_ptr);
  static void generate_transformer_U(const alps::params &parms,
                                     boost::shared_ptr<FourierTransformer> &fourier_ptr,
                                     const std::vector<double> &densities);
  static void generate_transformer_U(const alps::params &parms,
                                     boost::shared_ptr<FourierTransformer> &fourier_ptr,
                                     const std::vector<double> &densities,
                                     const std::vector<double> &magnetization);
  //void setmu(double mu){ mu_=mu;}
  
  double moments(const unsigned int z, const unsigned int i, const unsigned int j, 
                 const unsigned int mom) const 
  {
    switch (mom) {
      case 0 :
        return c1_[z][i][j];
      case 1 :
        return c2_[z][i][j];
      case 2 :
        return c3_[z][i][j];
      default :
        throw std::invalid_argument("no such coefficient stored in FourierTransformer");
    }
  }
  double SE_moments(const unsigned int z, const unsigned int i, const unsigned int j, 
                 const unsigned int mom) const 
  {
    switch (mom) {
      case 0 :
        return Sc0_[z][i][j];
      case 1 :
        return Sc1_[z][i][j];
      case 2 :
        return Sc2_[z][i][j];
      default :
        throw std::invalid_argument("no such coefficient stored in FourierTransformer");
    }
  }
  
  
  
  
  
protected:
  
  double beta_;
  //double mu_;
  std::vector<std::vector<std::vector<double> > > c1_;
  std::vector<std::vector<std::vector<double> > > c2_;
  std::vector<std::vector<std::vector<double> > > c3_;
  std::vector<std::vector<std::vector<double> > > Sc0_; //coefficients for the self-energy
  std::vector<std::vector<std::vector<double> > > Sc1_;  
  std::vector<std::vector<std::vector<double> > > Sc2_;  
  
  
};



class SimpleG0FourierTransformer : public FourierTransformer
{
public: 
  SimpleG0FourierTransformer(const double beta, const double mu, const double h, const int n_flavor, 
                             const std::vector<double>& eps, const std::vector<double>& epssq)
  : FourierTransformer(beta, n_flavor, 1)
  {
    for(int f=0; f<n_flavor; f++) {
      int s = f % 2 ? -1 : 1;
      double mub = mu + s*h;
      c1_[f][0][0] = 1.;
      c2_[f][0][0] = eps[f] - mub;
      c3_[f][0][0] = epssq[f] - 2*mub*eps[f] + mub*mub;	
    }
  }
};



class ClusterG0FourierTransformer : public FourierTransformer
{
public: 
  ClusterG0FourierTransformer(const double beta, const double mu, const double h, const int n_flavor, const int n_site, 
                              const std::vector<std::vector<double> >& eps, const std::vector<std::vector<double> >& epssq)
  : FourierTransformer(beta, n_flavor, n_site)
  { 
    if (n_site%2!=0 && h!=0) 
      boost::throw_exception(std::logic_error("FourierTransformer: cannot handle antiferromagnetism with odd cluster size\n"));
    for(int f=0; f<n_flavor; f++) {
      int s = f % 2 ? -1 : 1;
      for(int i=0; i<n_site; i++) {
        c1_[f][i][i] = 1.;
        c2_[f][i][i] = eps[f][i] - mu;
        c3_[f][i][i] = epssq[f][i] - 2*mu*eps[f][i] + mu*mu +h*h/4.;
      }
      if (n_site%2==0) {
        for(int i=0; i<n_site; i+=2) {
          c2_[f][i+1][i] = s*h/2;
          c2_[f][i][i+1] = c2_[f][i+1][i];
          c3_[f][i+1][i] = s*h/2*(eps[f][i] + eps[f][i+1] - 2*mu);
          c3_[f][i][i+1] = c3_[f][i+1][i];
        }
      }
    }
   /* for(int f=0;f<n_flavor;f++){
      for(int i=0;i<n_site;i++){
        for(int j=0;j<n_site;j++){
          std::cout<<f<<i<<j<<" "<<c1_[f][i][j]<<" "<<c2_[f][i][j]<<" "<<c3_[f][i][j]<<std::endl;
        }
      }
    }*/
  }
};



class GFourierTransformer : public FourierTransformer
{
public: 
  GFourierTransformer(const double beta, const double mu, const double U, const int n_flavor, const int n_site, 
                      const std::vector<double>& densities, 
                      const std::vector<std::vector<double> >& eps, const std::vector<std::vector<double> >& epssq)
  : FourierTransformer(beta, n_flavor, n_site)  
  {
    for(int f=0;f<n_flavor;++f){
      int fbar = f%2==0 ? f+1 : f-1;
      //std::cerr << "dens: " << densities[fbar] << std::endl;
      for(int i=0;i<n_site;++i){
        c1_[f][i][i] = 1.;
        c2_[f][i][i] = eps[f][i] - mu + U*densities[fbar]; 
        c3_[f][i][i] = epssq[f][i] - 2.*mu*eps[f][i] + mu*mu 
        + 2.*U*densities[fbar]*(eps[f][i]-mu) + U*U*densities[fbar];
        Sc0_[f][i][i] = U * (densities[fbar]-0.5);
        Sc1_[f][i][i] = U*U * densities[fbar] * (1-densities[fbar]);
        //std::cout << "eps: " << f << " " << i << " " << eps[f][i] << "\n";
        Sc2_[f][i][i] = 0;/*U*U * densities[fbar] * (1-densities[fbar])
                           * ( U * (1-densities[fbar]) + (eps[f][i]-mu) );*/
      }
    }
  }
};




class SimpleAntiFerroGFourierTransformer : public FourierTransformer
{
public: 
  SimpleAntiFerroGFourierTransformer(const double beta, const double /*mu*/, const double /*h*/, const double /*U*/, const int n_flavor, 
                                     const std::vector<double>& /*densities*/, const std::vector<double>& /*magnetization*/,  
                                     const std::vector<double>& /*eps*/, const std::vector<double>& /*epssq*/)
  : FourierTransformer(beta, n_flavor, 1)
  {
    std::cerr << "to be implemented\n";
    exit(1);
    for(int f=0;f<n_flavor;++f){
      //int fbar = f%2==0 ? f+1 : f-1;
      //int s = f % 2 ? -1 : 1;
      //double mub = mu + s*h;
      c1_[f][0][0] = 1.;
      c2_[f][0][0] = 0;
      c3_[f][0][0] = 0;
      Sc0_[f][0][0] = 0;
      Sc1_[f][0][0] = 0;
      Sc2_[f][0][0] = 0;
    }
  }
};



class ClusterAntiFerroGFourierTransformer : public FourierTransformer
{
public: 
  ClusterAntiFerroGFourierTransformer(const double beta, const double mu, const double h, const double U, 
                                      const int n_flavor, const int n_site, 
                                      const std::vector<double>& density, const std::vector<double>& magnetization, 
                                      const std::vector<std::vector<double> >& eps, const std::vector<std::vector<double> >& epssq)
  : FourierTransformer(beta, n_flavor, n_site)
  { 
    for(int f=0;f<n_flavor;++f){
      int fbar = f%2==0 ? f+1 : f-1;
      int s = f % 2 ? -1 : 1;
      for(int i=0;i<n_site;++i){
        c1_[f][i][i] = 1.;
        c2_[f][i][i] = eps[f][i] - mu + U*density[fbar]; 
        c3_[f][i][i] = epssq[f][i] - 2.*mu*eps[f][i] + mu*mu + h*h/4 
        + 2.*U*density[fbar]*(eps[f][i]-mu) + U*U*density[fbar] + U*s*h*magnetization[fbar];
        Sc0_[f][i][i] = U * (density[fbar]-0.5);
        Sc1_[f][i][i] = U*U * (density[fbar]*(1-density[fbar]) - magnetization[fbar]*magnetization[fbar]);
        /*double m = magnetization[fbar];
         double n = density[fbar];
         double e = eps[f][i];
         double eq = i%2 ? eps[f][i-1] : eps[f][i+1];*/
        Sc2_[f][i][i] = 0.;/*U*U*n*(1-n) * (U*(1-n) + (e-mu))
                            + U*U*m * (2*U*(n-1)*m + U*m*n - (eq-mu)*n + eq)
                            + U*s*h*m* ( U*(0.5-n) + (e-eq) );*/
      }
      for(int i=0; i<n_site; i+=2) {
        c2_[f][i+1][i] = s*h/2 + U*magnetization[fbar];
        c2_[f][i][i+1] = c2_[f][i+1][i];
        c3_[f][i+1][i] = (U*magnetization[fbar]+s*h/2)*(eps[f][i] + eps[f][i+1] - 2*mu)
        + U*U*magnetization[fbar] + U*s*h*density[fbar];
        c3_[f][i][i+1] = c3_[f][i+1][i];
        Sc0_[f][i+1][i] = U*magnetization[fbar];
        Sc0_[f][i][i+1] = Sc0_[f][i+1][i];
        Sc1_[f][i+1][i] = U*U*magnetization[fbar]*(1-2*density[fbar]);
        /*double m = magnetization[fbar];
         double n = density[fbar];
         double e = eps[f][i];
         double eq = i%2 ? eps[f][i-1] : eps[f][i+1];*/
        Sc2_[f][i][i+1] = 0;/*U*U*U * (m*m*m + m*(3*n*n-4*n+1)) 
                             - 0.5*U*U*s*h* (n*(1-n) + m*m)
                             +U*U*(((2*mu-eq-e)*n-mu)*m + 0*(e+eq));*/
        Sc2_[f][i][i+1] = Sc2_[f][i+1][i];	
      }
    }
    /*std::ofstream forstr("for.dat");
     for(int f=0;f<n_flavor;++f){
     std::cout << "density: " << density[f] << std::endl; 
     for(int s=0;s<n_site;s+=1){
     for (int p1=0; p1<=0; ++p1) {
     for (int p2=0; p2<=0; ++p2) {
     for(uint k=0; k<100; k++){
     std::complex<double> iw(0,(2*k+1)*M_PI/beta_);
     forstr << iw.imag() << "\t" << f << "\t" << s+p1 << "\t" << s+p2 << "\t"
     << f_omega(iw, c1_[f][s+p1][s+p2],c2_[f][s+p1][s+p2], c3_[f][s+p1][s+p2]).real() << "\t"
     << f_omega(iw, c1_[f][s+p1][s+p2],c2_[f][s+p1][s+p2], c3_[f][s+p1][s+p2]).imag()
     << std::endl;
     }
     }
     }
     }
     }*/
  }
  
  virtual void append_tail(matsubara_green_function_t& G_omega, const matsubara_green_function_t& G0_omega,
                           const unsigned int nfreq_measured) const{ throw std::logic_error("this is not implemented.");}
  
};




class FFunctionFourierTransformer:public FourierTransformer
{
public:
  FFunctionFourierTransformer(double beta, double mu, double epsilonsq_av, int n_flavor, int n_site)
  :FourierTransformer(beta, n_flavor, n_site){
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
  
  virtual void time_to_frequency_ft(const itime_green_function_t &G_tau, matsubara_green_function_t &G_omega) const{
    std::cout<<"implement 2nd derivatives for F function before using this!"<<std::endl;
    FourierTransformer::time_to_frequency_ft(G_tau, G_omega);
  }
  virtual void time_to_frequency_ft(const nambu_itime_green_function_t &G_tau, nambu_matsubara_green_function_t &G_omega) const{ throw std::logic_error("not implemented.");}
  
  virtual void frequency_to_time_ft(itime_green_function_t &G_tau, const matsubara_green_function_t &G_omega) const{
    FourierTransformer::frequency_to_time_ft(G_tau, G_omega);
    std::cout<<"correcting for F function"<<std::endl;
    for(unsigned int j=0;j<G_tau.nflavor();++j){
      G_tau(G_tau.ntime()-1,j)=-epsilonsq_av_-G_tau(0,j);
    }
  }
  virtual void frequency_to_time_ft(nambu_itime_green_function_t &G_tau, const nambu_matsubara_green_function_t &G_omega) const{throw std::logic_error("not implemented");}
  
  /*  virtual void frequency_to_time_ft(green_function<alps::RealObsevaluator> &G_tau, 
   const green_function<alps::RealObsevaluator> &G_omega_real,
   const green_function<alps::RealObsevaluator> &G_omega_imag) const{
   
   std::cout<<"backward ft: "<<&G_tau<<std::endl;
   std::cout<<"backward ft: "<<&G_omega_real<<std::endl;
   std::cout<<"backward ft: "<<&G_omega_imag<<std::endl;
   std::cout<<"not implemented for obsevaluators. exiting."<<std::endl;
   abort();
   } //not implemented.*/
  
  virtual ~FFunctionFourierTransformer(){}
  
private:
  double epsilonsq_av_;
};



class FittingFourierTransformer : public FourierTransformer
{
  //a simple Fourier Transformer which fits c1_/w to the tail of the matsubara green function
public:
  
  FittingFourierTransformer(double beta, int n_flavor, int n_site)
  : FourierTransformer(beta, n_flavor, n_site)
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
  
  virtual void time_to_frequency_ft(const itime_green_function_t &G_tau, matsubara_green_function_t &G_omega) 
  {
    std::cout<<"implement 2nd derivatives for F function before using this!"<<std::endl;
    FourierTransformer::time_to_frequency_ft(G_tau, G_omega);
  }
  virtual void time_to_frequency_ft(const nambu_itime_green_function_t &G_tau, nambu_matsubara_green_function_t &G_omega) const{throw std::logic_error("not implemented");}
  
  virtual void frequency_to_time_ft(itime_green_function_t &G_tau, const matsubara_green_function_t &G_omega) const;
  virtual void frequency_to_time_ft(nambu_itime_green_function_t &G_tau, const nambu_matsubara_green_function_t &G_omega) const{throw std::logic_error("not implemented");}
  
  /*  virtual void frequency_to_time_ft(green_function<alps::RealObsevaluator> &G_tau, 
   const green_function<alps::RealObsevaluator> &G_omega_real,
   const green_function<alps::RealObsevaluator> &G_omega_imag) const
   {
   
   std::cout<<"backward ft: "<<&G_tau<<std::endl;
   std::cout<<"backward ft: "<<&G_omega_real<<std::endl;
   std::cout<<"backward ft: "<<&G_omega_imag<<std::endl;
   std::cout<<"not implemented for obsevaluators. exiting."<<std::endl;
   abort();
   } //not implemented.*/
  
  virtual ~FittingFourierTransformer(){}
  
};




#endif
