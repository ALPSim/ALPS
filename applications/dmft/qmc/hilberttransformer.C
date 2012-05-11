/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2007 by Philipp Werner <werner@comp-phys.org>,
 *                         Matthias Troyer <troyer@comp-phys.org>,
 *                         Emanuel Gull <gullc@comp-phys.org>,
 *                         Sebastian Fuchs <fuchs@comp-phys.org>
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

/* $Id: hilberttransformer.C 338 2009-01-20 01:22:05Z fuchs $ */



#include "hilberttransformer.h"
#include <functional>
#include <math.h>

itime_green_function_t HilbertTransformer::symmetrize(const itime_green_function_t& G_tau, const bool symmetrization) const
{
  itime_green_function_t G(G_tau);
  if (symmetrization) {
    assert(G.nflavor()%2==0);
    for(spin_t flavor=0;flavor<G.nflavor(); flavor+=2){
      for(itime_index_t tau=0;tau<G.ntime();++tau){
        G(tau, flavor  )=0.5*(G(tau, flavor)+G(tau, flavor+1));
        G(tau, flavor+1)=G(tau, flavor);
      }
    }
  }
  return G;
}


itime_green_function_t HilbertTransformer::initial_G0(const alps::Parameters& parms) 
{
  throw std::logic_error("not implemented - specify your Hilbert transformer");
}




itime_green_function_t SemicircleHilbertTransformer::operator()(const itime_green_function_t& G_tau, 
                                                                double mu, double h, double beta)
{
  double tsq=t_*t_;
  matsubara_green_function_t G_omega(G_tau.ntime()-1, G_tau.nsite(), G_tau.nflavor());
  matsubara_green_function_t G0_omega(G_tau.ntime()-1, G_tau.nsite(), G_tau.nflavor());
  itime_green_function_t G0_tau(G_tau.ntime(), G_tau.nsite(), G_tau.nflavor());
  std::cerr << "SemicircleHilbertTransformer::operator(): Fouriertransformation of G_tau at this point is requested. Densities needed!\n";
  exit(1);
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  //FourierTransformer::generate_transformer(parms, fourier_ptr); ????
  fourier_ptr->forward_ft(G_tau, G_omega);
  //std::cout<<"G omega real: "<<std::endl;
  print_real_green_matsubara(std::cout, G_omega, beta);
  for(unsigned i=0; i<G_omega.nfreq(); i++) {
    std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
    for(spin_t flavor=0;flavor<G_omega.nflavor(); flavor++){
      std::complex<double> zeta = iw + mu + (flavor%2 ? -h : h);
      G0_omega(i, flavor) = 1./(zeta - tsq*G_omega(i, flavor%2==0?flavor+1:flavor-1));
    }
  }
  //std::cout<<"symmetrized G0 omega real: "<<std::endl;
  print_real_green_matsubara(std::cout, G0_omega, beta);
  fourier_ptr->backward_ft(G0_tau, G0_omega);
  return G0_tau;
}




itime_green_function_t SemicircleHilbertTransformer::initial_G0(const alps::Parameters& parms) 
{
  std::cout<<"calculating initial G0_tau for semicircular DOS"<<std::endl;
  int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  int n_time=boost::lexical_cast<int>(parms["N"]);
  int n_flavor=parms.value_or_default("FLAVORS", 2);
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  double h = static_cast<double>(parms.value_or_default("H_INIT",0.));
  std::vector<double> tsq(n_flavor, t_*t_);
  for(int i=0;i<n_flavor;++i){
    std::stringstream ti; ti<<"t"<<i/2;
    if(parms.defined(ti.str())){
      double t=static_cast<double>(parms[ti.str()]);
      tsq[i]=t*t;
      std::cout<<"for flavor: "<<i<<" using bw: "<<t<<std::endl;
    }
  }
  matsubara_green_function_t G0_omega(n_matsubara,n_flavor);  
  itime_green_function_t G0_tau(n_time+1, n_flavor);
  for(spin_t flavor=0;flavor<G0_omega.nflavor(); flavor++) {
    for(unsigned i=0; i<G0_omega.nfreq(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      std::complex<double> zeta = iw+mu+(flavor%2 ? -h : h);
      G0_omega(i, flavor) = (zeta - sqrt(zeta*zeta-4*tsq[flavor]))/(2*tsq[flavor]);
    }
  }
  G0_omega.write("G0_omega_input");
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  fourier_ptr->backward_ft(G0_tau, G0_omega);
  G0_tau.write("G0_tau_input");
  return G0_tau;
}




matsubara_green_function_t FrequencySpaceHilbertTransformer::initial_G0(const alps::Parameters& parms) 
{
  std::cout<<"reading initial G0_omega"<<std::endl;
  int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  int n_orbital=parms.value_or_default("FLAVORS", 2);
  matsubara_green_function_t G0_omega(n_matsubara, n_orbital);
  if (parms.defined("INSULATING")) {
    std::cout<<"calculating initial G0_omega"<<std::endl;
    int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
    int n_time=boost::lexical_cast<int>(parms["N"]);
    int n_flavor=parms.value_or_default("FLAVORS", 2);
    double beta = static_cast<double>(parms["BETA"]);
    double mu = static_cast<double>(parms["MU"]);
    double h = static_cast<double>(parms["H_INIT"]);
    matsubara_green_function_t G0_omega(n_matsubara,n_flavor);
    for(unsigned i=0; i<G0_omega.nfreq(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      for(spin_t flavor=0;flavor<G0_omega.nflavor(); flavor++) {
        std::complex<double> zeta = iw+mu+(flavor%2 ? -h : h);
        G0_omega(i, flavor) = 1./zeta;
      }
    }
    G0_omega.write("G0_omega_input");
    boost::shared_ptr<FourierTransformer> fourier_ptr;
    FourierTransformer::generate_transformer(parms, fourier_ptr);
    itime_green_function_t G0_tau(n_time,n_flavor);
    fourier_ptr->backward_ft(G0_tau, G0_omega);
    G0_tau.write("G0_tau_input");
    return G0_omega;
  }
  else {
    std::ifstream check(parms["G0OMEGA_INPUT"].c_str());
    int n_time=boost::lexical_cast<int>(parms["N"]);
    int n_flavor=parms.value_or_default("FLAVORS", 2);
    if(!check.is_open()){
      std::cout<<"could not open inital G0 file "<<parms["G0OMEGA_INPUT"]<<"; proceeding with free solution. "<<std::endl;
      matsubara_green_function_t G_omega(n_matsubara, n_orbital);
      for(int i=0;i<n_matsubara;++i){
        for(int f=0;f<n_orbital;++f){
          G_omega(i,f)=1;
          G0_omega(i,f)=1;
        }
      }
      double h     =parms.value_or_default("H_INIT",0.);
      double mu    =parms.value_or_default("MU",0.);
      double beta  =parms.value_or_default("BETA",0.);
      G0_omega=this->operator()(G_omega, G0_omega, mu, h, beta);
      G0_omega.write("G0_omega_input");
      boost::shared_ptr<FourierTransformer> fourier_ptr;
      FourierTransformer::generate_transformer(parms, fourier_ptr);
      itime_green_function_t G0_tau(n_time,n_flavor);
      fourier_ptr->backward_ft(G0_tau, G0_omega);
      G0_tau.write("G0_tau_input");
    }
    else
      G0_omega.read(parms["G0OMEGA_INPUT"].c_str());
  }
  return G0_omega;
}




matsubara_green_function_t FSSemicircleHilbertTransformer::operator()(const matsubara_green_function_t& G_omega, 
                                                                      matsubara_green_function_t &G0_omega_ignored, 
                                                                      const double mu, const double h, const double beta)
{
  std::cout<<"Semicircular Hilbert Transform: using: mu:"<<mu<<" h: "<<h<<"beta: "<<beta<<" "<<G0_omega_ignored.nfreq()<<std::endl;
  double tsq=t_*t_;
  matsubara_green_function_t G0_omega(G_omega);
  //formula according to review, p. 61, formula 221. 
  if(G_omega.nflavor()==1){ //special case. 
    for(unsigned i=0; i<G_omega.ntime(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      G0_omega(i,0) =1./(iw + mu - tsq*G_omega(i,0));
    }
  }
  else{
    assert(G_omega.nflavor()%2==0);
    for(unsigned flavor=0;flavor<G_omega.nflavor();flavor+=2){
      for(unsigned i=0; i<G_omega.nfreq(); i++) {
        std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
        G0_omega(i,flavor  ) =1./(iw + mu -h - tsq*G_omega(i,flavor+1)); 
        G0_omega(i,flavor+1) =1./(iw + mu +h - tsq*G_omega(i,flavor  )); 
      }
    }
  }
  return G0_omega;
}




matsubara_green_function_t FSSemicircleHilbertTransformer::initial_G0(const alps::Parameters& parms)
{
  if (parms.defined("G0OMEGA_INPUT") || parms.defined("INSULATING"))
    return FrequencySpaceHilbertTransformer::initial_G0(parms);
  else {
    std::cout<<"calculating initial G0_omega"<<std::endl;
    int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
    int n_time=boost::lexical_cast<int>(parms["N"]);
    int n_flavor=parms.value_or_default("FLAVORS", 2);
    double beta = static_cast<double>(parms["BETA"]);
    double mu = static_cast<double>(parms["MU"]);
    double h = static_cast<double>(parms["H_INIT"]);
    double tsq=t_*t_;
    matsubara_green_function_t G0_omega(n_matsubara,n_flavor);
    for(unsigned i=0; i<G0_omega.nfreq(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      for(spin_t flavor=0;flavor<G0_omega.nflavor(); flavor++) {
        std::complex<double> zeta = iw+mu+(flavor%2 ? -h : h);
        std::complex<double> tmp = sqrt(zeta*zeta-4*tsq);
        tmp *= tmp.imag()<0 ? -1 : 1;
        G0_omega(i, flavor) = (zeta - tmp)/(2*tsq);
        //std::cout << zeta << " " << sqrt(zeta*zeta) << std::endl;
      }
    }
    G0_omega.write("G0_omega_input");
    boost::shared_ptr<FourierTransformer> fourier_ptr;
    FourierTransformer::generate_transformer(parms, fourier_ptr);
    itime_green_function_t G0_tau(n_time,n_flavor);
    fourier_ptr->backward_ft(G0_tau, G0_omega);
    G0_tau.write("G0_tau_input");
    return G0_omega;
  }
}




FSDOSHilbertTransformer::FSDOSHilbertTransformer(alps::Parameters& params)
{
  if(!params.defined("DOSFILE")) throw std::invalid_argument("please define DOSFILE parameter!");
  std::ifstream dos_file(params["DOSFILE"].c_str());
  if(!dos_file.good()) throw std::runtime_error("DOSFILE is not good!");
  double eps, d;
  std::cout<<"using density of states "<<params["DOSFILE"]<<std::endl;
  while(dos_file>>eps>>d){
    epsilon.push_back(eps);
    dos.push_back(d);
    //std::cout<<eps<<" "<<d<<std::endl;
  }
  if(dos.size()%2!=1){ throw std::runtime_error("please use odd number of DOS points (due to Simpson integration)"); }
  //normalize DOS to one:
  double s=0;
  for(unsigned i=1;i<epsilon.size()-2;i+=2){
    s+=4.*dos[i]+2*dos[i+1];                   // Assuming equidistant epsilon points
  }
  s+=dos[0]+dos[epsilon.size()-1]+4.*dos[epsilon.size()-2];
  s/=(3.);
  std::cout<<"normalization constant: "<<s<<std::endl;
  std::cout<<"step size h is: " << (epsilon[1]-epsilon[0]) << "   (Note: Assuming equidistant energy intervals.)" << std::endl;
  for(unsigned i=0;i<epsilon.size();++i){
    dos[i]/=s;
  }
  double S=0;
  for(unsigned i=0;i<epsilon.size()-1;++i){
    S+=dos[i];
  }
  std::cout<<"check: total sum is: "<<S<<std::endl;
  //find second moment of band structure:
  double epssqav=0;
  for(unsigned i=1;i<epsilon.size()-2;i+=2){
    epssqav+=4.*(epsilon[i]*epsilon[i]*dos[i])+2*(dos[i+1]*epsilon[i+1]*epsilon[i+1]);
  } 
  epssqav+=dos[0]*epsilon[0]*epsilon[0]+dos[epsilon.size()-1]*epsilon[epsilon.size()-1]*epsilon[epsilon.size()-1]+4.*dos[epsilon.size()-2]*epsilon[epsilon.size()-2]*epsilon[epsilon.size()-2];
  epssqav/=(3.);
  std::cout<<"second moment of band structure: "<<epssqav<<std::endl;
  
  // set parameters for the FourierTransformer if they are not set by the user
  bool is_set = params.defined("EPSSQAV");
  unsigned int n_flavors = params.value_or_default("FLAVORS", 2);
  for (unsigned int i=0; i < n_flavors; i++) {
    if (params.defined("EPS_"+boost::lexical_cast<std::string>(i)) || params.defined("EPSSQ_"+boost::lexical_cast<std::string>(i)))
      is_set = true;
  }
  if (! is_set) {
    std::cout << "FSDOSHilbertTransformer: writing the values for bandstructure-parameters needed in FourierTransformer." << std::endl;
    params["EPSSQAV"] = epssqav;
    for (unsigned int i=0; i < n_flavors; i++) {
      params["EPS_"+boost::lexical_cast<std::string>(i)] = 0;
      params["EPSSQ_"+boost::lexical_cast<std::string>(i)] = epssqav;
    }
  }
}




matsubara_green_function_t FSDOSHilbertTransformer::operator()(const matsubara_green_function_t &G_omega, 
                                                               matsubara_green_function_t &G0_omega, 
                                                               const double mu, const double h, const double beta)
{
  assert(G_omega.nsite()==1);
  if(h!=0){throw std::invalid_argument(" we don't yet handle magnetic fields!"); }
  
  matsubara_green_function_t G_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
  matsubara_green_function_t G0_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
  matsubara_green_function_t sigma(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
  for(spin_t f=0;f<G_omega.nflavor();++f){
    for(frequency_t w=0;w<G_omega.nfreq();++w){
      std::complex<double> g=0;
      sigma(w,f)= 1./G0_omega(w,f)-1./G_omega(w,f); 
      std::complex<double> A=std::complex<double>(0, (2*w+1)*M_PI/beta) + mu -sigma(w, f);
      for(unsigned n=0;n<dos.size()-1;++n){
        g+=dos[n]/(A-epsilon[n]); //go back to higher order integrator!!
      }
      G_omega_new(w, f)=g;
    }
  }
  if(G_omega.nflavor() !=2){throw std::logic_error("don't know how to handle != 2 flavors.");}
  for(spin_t f=0;f<G_omega.nflavor();++f){
    for(frequency_t w=0;w<G_omega.nfreq();++w){
      G0_omega_new(w,f)=1./(sigma(w,f)+1./G_omega_new(w,f));
    }
  }
  return G0_omega_new;
}




matsubara_green_function_t AFM_FSDOSHilbertTransformer::operator()(const matsubara_green_function_t &G_omega,
                                                                   matsubara_green_function_t &G0_omega, 
                                                                   const double mu, const double h, const double beta)
{
  if(G_omega.nsite()!=1){ throw std::invalid_argument("don't know how to handle n_site >1. Use Sebastian's cluster loop for that. aborting.");}
  if(G_omega.nflavor()!=2){ throw std::invalid_argument("don't know how to handle n_flavor!=2. aborting.");}
  std::cout<<"AFM FS DOS Hilbert Transformer"<<std::endl;
  matsubara_green_function_t Sigma(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
  matsubara_green_function_t G_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
  matsubara_green_function_t G0_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
  //compute self energy
  for(frequency_t w=0;w<G_omega.nfreq();++w){
    double wn=(2.*w+1)*M_PI/beta;
    for(spin_t f=0;f<G_omega.nflavor();++f){
      Sigma(w,f)=1./G0_omega(w,f)-1./G_omega(w,f);
    }
    //std::cout<<"sigma: "<<Sigma(w,0)<<" "<<Sigma(w,1)<<std::endl;
    std::complex<double>zeta_A=std::complex<double>(mu-h,wn)-Sigma(w,0);
    std::complex<double>zeta_B=std::complex<double>(mu+h,wn)-Sigma(w,1);
    
    //Simpson: integrate dos(e)/(zeta_A zeta_B-e^2)
    if(dos.size()%2 !=1){throw std::runtime_error("for Simpson precision: use a DOS with odd number of integration points. "); }
    std::complex<double> I=0;
    for(unsigned i=1;i<dos.size()-2;i+=2){
      I+=4.*integrand(i,dos,epsilon,zeta_A,zeta_B)+2.*integrand(i+1,dos,epsilon,zeta_A,zeta_B);
    }
    I+=integrand(0,dos,epsilon,zeta_A,zeta_B)+integrand(dos.size()-1,dos,epsilon,zeta_A,zeta_B)+4.*integrand(dos.size()-2,dos,epsilon,zeta_A,zeta_B);
    I/=3.;
    //compute the new G's:
    G_omega_new(w,0)=zeta_B*I; //formula 97 in review
    G_omega_new(w,1)=zeta_A*I;
    G0_omega_new(w,0)=1./(1./G_omega_new(w,0)+Sigma(w,0));
    G0_omega_new(w,1)=1./(1./G_omega_new(w,1)+Sigma(w,1));
  }
  return G0_omega_new;
}






matsubara_green_function_t TwoDAFMHilbertTransformer::operator()(const matsubara_green_function_t & G_omega, 
                                                                 matsubara_green_function_t &G0_omega, 
                                                                 const double mu, const double h, const double beta)
{
  //compute sigma 
  matsubara_green_function_t sigma(G_omega);
  matsubara_green_function_t G(G_omega);
  for(unsigned i=0;i<G_omega.nfreq();++i){
    sigma(i,0)=1./G0_omega(i,0)-1./G_omega(i,0);
    sigma(i,1)=1./G0_omega(i,1)-1./G_omega(i,1);
  }
  //perform integration over bz
  for(unsigned i=0;i<G_omega.nfreq();++i){
    G(i,0)=0.;
    G(i,1)=0.;
    std::complex<double> iomegan(0,(2*i+1)*M_PI/beta);
    std::complex<double> zeta_A=iomegan + mu -h -sigma(i,0);
    std::complex<double> zeta_B=iomegan + mu +h -sigma(i,1);
    afm_integrand integrand(t_,tprime_,zeta_A, zeta_B);
    std::complex<double> I=twodsimpson(integrand, -M_PI, -M_PI, M_PI, M_PI, L_);
    G(i,0)=1./(4*M_PI*M_PI)*zeta_B*I;
    G(i,1)=1./(4*M_PI*M_PI)*zeta_A*I;
  }
  //compute the new G0 
  for(unsigned i=0;i<G_omega.nfreq();++i){
    G0_omega(i,0)=1./(sigma(i,0)+1./G(i,0));
    G0_omega(i,1)=1./(sigma(i,1)+1./G(i,1));
  }
  return G0_omega;
}




