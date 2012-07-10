/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2007 by Philipp Werner <werner@comp-phys.org>,
 *                            Matthias Troyer <troyer@comp-phys.org>,
 *                            Emanuel Gull <gullc@comp-phys.org>,
 *                            Sebastian Fuchs <fuchs@comp-phys.org>
 *               2012      by Jakub Imriska <jimriska@phys.ethz.ch>
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
      G0_omega(i, flavor) = 1./(zeta - bethe_parms.tsq(flavor)*G_omega(i, flavor%2==0?flavor+1:flavor-1));
    }
  }
  //std::cout<<"symmetrized G0 omega real: "<<std::endl;
  print_real_green_matsubara(std::cout, G0_omega, beta);
  fourier_ptr->backward_ft(G0_tau, G0_omega);
  return G0_tau;
}




itime_green_function_t SemicircleHilbertTransformer::initial_G0(const alps::Parameters& parms) 
{
  std::cout<<"SemicircularHilbertTransformer::initial_G0: ";
  int n_time=boost::lexical_cast<int>(parms["N"]);
  int n_flavor=parms.value_or_default("FLAVORS", 2);
  itime_green_function_t G0_tau(n_time+1, n_flavor);
  
  if (parms.defined("G0TAU_INPUT") && parms["G0TAU_INPUT"].length()>0) {
	std::cout<<"reading initial G0_tau"<<std::endl;
    std::ifstream check(parms["G0TAU_INPUT"].c_str());
    if(!check.is_open()) {
      std::cerr << "ERROR: could not open inital G0 file "<<parms["G0TAU_INPUT"]<<std::endl;
      throw std::runtime_error("SemicircleHilbertTransformer::initial_G0: could not open inital G0 file");
    }
    else
      G0_tau.read(parms["G0TAU_INPUT"].c_str());
  }
  else {
    FSSemicircleHilbertTransformer hilbert(parms);
    boost::shared_ptr<FourierTransformer> fourier_ptr;
    FourierTransformer::generate_transformer(parms, fourier_ptr,true);   // with H_INIT (if provided)
    fourier_ptr->backward_ft(G0_tau, hilbert.initial_G0(parms));
  }
  
  if (parms.defined("G0TAU_input"))
    G0_tau.write((parms["G0TAU_input"]).c_str());
  
  return G0_tau;
}




matsubara_green_function_t FrequencySpaceHilbertTransformer::initial_G0(const alps::Parameters& parms) 
{
  std::cout<<"FrequencySpaceHilbertTransformer::initial_G0: ";
  int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  int n_time=boost::lexical_cast<int>(parms["N"]);
  int n_orbital=parms.value_or_default("FLAVORS", 2);
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  double h = (parms.defined("H_INIT") ? static_cast<double>(parms["H_INIT"]) : static_cast<double>(parms.value_or_default("H",0.)));
  matsubara_green_function_t G0_omega(n_matsubara, n_orbital);

  if (parms.defined("INSULATING")) {
    std::cout<<"calculating insulating initial G0_omega"<<std::endl;
    for(unsigned i=0; i<G0_omega.nfreq(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      for(spin_t flavor=0;flavor<G0_omega.nflavor(); flavor++) {
        std::complex<double> zeta = iw+mu+(flavor%2 ? -h : h);
        G0_omega(i, flavor) = 1./zeta;
      }
    }
  }
  else if (parms.defined("G0OMEGA_INPUT") && parms["G0OMEGA_INPUT"].length()>0) {
	std::cout<<"reading initial G0_omega"<<std::endl;
    std::ifstream check(parms["G0OMEGA_INPUT"].c_str());
    if(!check.is_open()) {
      std::cerr << "ERROR: could not open inital G0 file "<<parms["G0OMEGA_INPUT"]<<std::endl;
      throw std::runtime_error("FrequencySpaceHilbertTransformer::initial_G0: could not open inital G0 file");
    }
    else
      G0_omega.read(parms["G0OMEGA_INPUT"].c_str());
  }
  else {
	std::cout<<"calculating non-interacting initial G0_omega"<<std::endl;
	if (parms.defined("DOSFILE") || parms.defined("TWODBS")) {
      matsubara_green_function_t G_omega(n_matsubara, n_orbital);
      for(int i=0;i<n_matsubara;++i){
        for(int f=0;f<n_orbital;++f){
          G_omega(i,f)=1;
          G0_omega(i,f)=1;
        }
      }
      G0_omega=this->operator()(G_omega, G0_omega, mu, h, beta);
    }
    else {   // Bethe lattice 
      BetheBandstructure bethe_parms(parms);
      for(unsigned i=0; i<n_matsubara; i++) {
        std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
        for(spin_t flavor=0;flavor<n_orbital; flavor++) {
          std::complex<double> zeta = iw+mu+(flavor%2 ? h : -h);
          std::complex<double> tmp = sqrt(zeta*zeta-4*bethe_parms.tsq(flavor));
          tmp *= tmp.imag()<0 ? -1 : 1;
          G0_omega(i, flavor) = (zeta - tmp)/(2*bethe_parms.tsq(flavor));
        }
      }
    }
  }
  
  if (parms.defined("G0OMEGA_input"))   // it is not needed to store it by default, as it will be stored in the 1st iteration as G0_omega_1, G0_omegareal_1
    G0_omega.write((parms["G0OMEGA_input"]).c_str());

  return G0_omega;
}




matsubara_green_function_t FSSemicircleHilbertTransformer::operator()(const matsubara_green_function_t& G_omega, 
                                                                      matsubara_green_function_t &G0_omega_ignored, 
                                                                      const double mu, const double h, const double beta)
{
  std::cout<<"Semicircular Hilbert Transform: using: mu:"<<mu<<" h: "<<h<<"beta: "<<beta<<" "<<G0_omega_ignored.nfreq()<<std::endl;
  matsubara_green_function_t G0_omega(G_omega);
  //formula according to review, p. 61, formula 221. 
  if(G_omega.nflavor()==1){ //special case. 
    for(unsigned i=0; i<G_omega.ntime(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      G0_omega(i,0) =1./(iw + mu - bethe_parms.tsq(0)*G_omega(i,0));
    }
  }
  else{
    assert(G_omega.nflavor()%2==0);
    for(unsigned flavor=0;flavor<G_omega.nflavor();flavor+=2){
      for(unsigned i=0; i<G_omega.nfreq(); i++) {
        std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
        G0_omega(i,flavor  ) =1./(iw + mu -h - bethe_parms.tsq(flavor)*G_omega(i,flavor+1)); 
        G0_omega(i,flavor+1) =1./(iw + mu +h - bethe_parms.tsq(flavor)*G_omega(i,flavor  )); 
      }
    }
  }
  return G0_omega;
}


void FrequencySpaceHilbertTransformer::SetBandstructureParms(alps::Parameters& parms, std::vector<double> eps, std::vector<double> epssq) {
  // set parameters for the FourierTransformer if none of them is set by the user
  bool is_set = parms.defined("EPSSQAV");
  unsigned int n_flavors = parms.value_or_default("FLAVORS", 2);
  for (unsigned int i=0; i < n_flavors; i++) {
    if (parms.defined("EPS_"+boost::lexical_cast<std::string>(i)) || parms.defined("EPSSQ_"+boost::lexical_cast<std::string>(i)))
      is_set = true;
  }
  if (! is_set) { //What is going on here?
    std::cout << "initialization of a the Hilbert transformer: writing the values for bandstructure-parameters needed in FourierTransformer." << std::endl;
    parms["EPSSQAV"] = epssq[0];
    for (unsigned int i=0; i < n_flavors; i++) {
      parms["EPS_"+boost::lexical_cast<std::string>(i)] = 0;
      parms["EPSSQ_"+boost::lexical_cast<std::string>(i)] = epssq[0];
    }
  }
}

FSDOSHilbertTransformer::FSDOSHilbertTransformer(alps::Parameters& parms)
{
  if(!parms.defined("DOSFILE")) throw std::invalid_argument("please define DOSFILE parameter!");
  std::ifstream dos_file(parms["DOSFILE"].c_str());
  if(!dos_file.good()) throw std::runtime_error("DOSFILE is not good!");
  unsigned int n_flavors = parms.value_or_default("FLAVORS", 2);
  if(n_flavors%2==1)n_flavors++;
  double eps, d;
  std::cout<<"using density of states "<<parms["DOSFILE"]<<std::endl;
  std::cout<<"(Note: Assuming equidistant energy intervals.)"<<std::endl;
  std::vector<double> dummyvec;
  for(int i = 0; i<n_flavors/2; i++){
    epsilon.push_back(dummyvec);
    dos.push_back(dummyvec);
  }
  while(dos_file>>eps>>d){
    epsilon[0].push_back(eps);
    dos[0].push_back(d);
    for(int i = 1; i<n_flavors/2; i++){
      if(!(dos_file>>eps>>d))throw std::runtime_error("DOSFILE is not good!");
      epsilon[i].push_back(eps);
      dos[i].push_back(d);
    }
  }

  double s;
  double S;
  std::vector<double> eps_;
  std::vector<double> epssq;

  for(int r = 1; r<n_flavors/2; r++){
    DOS_integrand integrand0(0, dos[r], epsilon[r]);
    s=(simpson_integrate(integrand0)).real();
    S=0;
    for(unsigned i=0;i<dos.size();++i){
      dos[r][i]/=s;
      S+=dos[r][i];
    }
    std::cout<<"check: total sum after normalization is: "<<S<<" (should close to 1)"<<std::endl;

    //find second moment of band structure: simpson integrate: dos[n]*epsilon[n]*epsilon[n]
    DOS_integrand integrand1(1,dos[r],epsilon[r]);
    DOS_integrand integrand2(2,dos[r],epsilon[r]);
    eps_.push_back((simpson_integrate(integrand1)).real());
    epssq.push_back((simpson_integrate(integrand2)).real());
    std::cout<<"Band "<<r<<": first moment of bandstructure: "<<eps_[r]<<std::endl;
    std::cout<<"Band "<<r<<": second moment of bandstructure: "<<epssq[r]<<std::endl;
  }
  
  SetBandstructureParms(parms,eps_,epssq);
}




matsubara_green_function_t FSDOSHilbertTransformer::operator()(const matsubara_green_function_t &G_omega, 
                                                               matsubara_green_function_t &G0_omega, 
                                                               const double mu, const double h, const double beta)
{
  if(G_omega.nsite()!=1){throw std::logic_error("FSDOSHilbertTransformer::operator(): don't know how to handle != 1 site.");}
  if(G_omega.nflavor()!=2){throw std::logic_error("FSDOSHilbertTransformer::operator(): don't know how to handle != 2 flavors.");}
  if(h!=0){throw std::invalid_argument("FSDOSHilbertTransformer::operator(): don't yet handle magnetic fields!");}
  std::cout<<"PM FS DOS Hilbert Transformer"<<std::endl;
  matsubara_green_function_t G_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
  matsubara_green_function_t Sigma(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
  
  for(spin_t f=0;f<G_omega.nflavor();++f){
    for(frequency_t w=0;w<G_omega.nfreq();++w){
      Sigma(w,f)= 1./G0_omega(w,f)-1./G_omega(w,f); 
      std::complex<double> A=std::complex<double>(0, (2*w+1)*M_PI/beta) + mu - Sigma(w,f);
      // simpson integrate: dos[n]/(A - epsilon[n])
      DOS_integrand integrand(dos[f],epsilon[f],A);
      G_omega_new(w, f)=simpson_integrate(integrand);
    }
  }
  
  for(spin_t f=0;f<G_omega.nflavor();++f){
    for(frequency_t w=0;w<G_omega.nfreq();++w){
      G0_omega(w,f)=1./(Sigma(w,f)+1./G_omega_new(w,f));
    }
  }
  return G0_omega;
}




matsubara_green_function_t AFM_FSDOSHilbertTransformer::operator()(const matsubara_green_function_t &G_omega,
                                                                   matsubara_green_function_t &G0_omega, 
                                                                   const double mu, const double h, const double beta)
{
  if(G_omega.nsite()!=1){ throw std::invalid_argument("AFM_FSDOSHilbertTransformer::operator(): don't know how to handle n_site >1. Use Sebastian's cluster loop for that. aborting.");}
  if(G_omega.nflavor()!=2){ throw std::invalid_argument("AFM_FSDOSHilbertTransformer::operator(): don't know how to handle n_flavor!=2. aborting.");}
  std::cout<<"AFM FS DOS Hilbert Transformer"<<std::endl;
  matsubara_green_function_t Sigma(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
  matsubara_green_function_t G_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());

  //compute self energy  *changed to many bands on 10th July, 2012 by René John Kerkdyk*
  for(frequency_t w=0;w<G_omega.nfreq();++w){
    double wn=(2.*w+1)*M_PI/beta;
    for(spin_t f=0;f<G_omega.nflavor();++f){
      Sigma(w,f)=1./G0_omega(w,f)-1./G_omega(w,f);
      if((f%2)==1){
	std::complex<double>zeta_A=std::complex<double>(mu-h,wn)-Sigma(w,f);
	std::complex<double>zeta_B=std::complex<double>(mu+h,wn)-Sigma(w,f);
	DOS_integrand integrand(dos[f],epsilon[f],zeta_A,zeta_B);
	std::complex<double> I = simpson_integrate(integrand);
	G_omega_new(w,f-1)=zeta_B*I; //formula 97 in review
	G_omega_new(w,f)=zeta_A*I;
	G0_omega(w,f-1)=1./(1./G_omega_new(w,f-1)+Sigma(w,f-1));
	G0_omega(w,f)=1./(1./G_omega_new(w,f)+Sigma(w,f));
      }
    }
    /*
    //std::cout<<"sigma: "<<Sigma(w,0)<<" "<<Sigma(w,1)<<std::endl; //OLD
    std::complex<double>zeta_A=std::complex<double>(mu-h,wn)-Sigma(w,0);
    std::complex<double>zeta_B=std::complex<double>(mu+h,wn)-Sigma(w,1);
    
    //Simpson: integrate dos(e)/(zeta_A zeta_B-e^2)
    DOS_integrand integrand(dos,epsilon,zeta_A,zeta_B);
    std::complex<double> I = simpson_integrate(integrand);

    //compute the new G's:
    G_omega_new(w,0)=zeta_B*I; //formula 97 in review
    G_omega_new(w,1)=zeta_A*I;
    G0_omega(w,0)=1./(1./G_omega_new(w,0)+Sigma(w,0));
    G0_omega(w,1)=1./(1./G_omega_new(w,1)+Sigma(w,1));
    */
  }
  return G0_omega;
}


TwoDHilbertTransformer::TwoDHilbertTransformer(alps::Parameters& parms) 
    : bandstruct_(parms),
      L_(parms.value_or_default("L",200))
{
  // set parameters for the FourierTransformer if they are not set by the user
  double epssq = bandstruct_.second_moment();
  std::cout<<"second moment of bandstructure: "<<epssq<<std::endl;
  std::vector<double> dummyvec(1,0.);
  std::vector<double> dummyvec2(1,epssq);
  if (epssq>0) SetBandstructureParms(parms,dummyvec,dummyvec2);
}


matsubara_green_function_t TwoDHilbertTransformer::operator() (const matsubara_green_function_t & G_omega, 
                                                                 matsubara_green_function_t &G0_omega, 
                                                                 const double mu, const double h, const double beta)
{
  if(G_omega.nsite()!=1){throw std::logic_error("TwoDHilbertTransformer::operator(): don't know how to handle != 1 site.");}
  if(G_omega.nflavor() !=2){throw std::logic_error("TwoDHilbertTransformer::operator(): don't know how to handle != 2 flavors.");}
  if(h!=0){throw std::invalid_argument("TwoDHilbertTransformer::operator(): don't yet handle magnetic fields!"); }
  std::cout<<"PM TwoD BS FS Hilbert Transformer"<<std::endl;
  matsubara_green_function_t Sigma(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
  matsubara_green_function_t G_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());

  //perform integration over bz
  for(unsigned f=0;f<G_omega.nflavor();++f){
    for(unsigned i=0;i<G_omega.nfreq();++i){
      Sigma(i,f)=1./G0_omega(i,f)-1./G_omega(i,f);
      std::complex<double> iomegan(0,(2*i+1)*M_PI/beta);
      std::complex<double> zeta_A=iomegan + mu - Sigma(i,f);
      std::complex<double> I;
      if (bandstruct_.single_band()) {
        pm_band_integrand integrand(bandstruct_,zeta_A);
        I=twodsimpson(integrand, bandstruct_.kx_min(), bandstruct_.ky_min(), bandstruct_.kx_max(), bandstruct_.ky_max(), L_) / bandstruct_.k_area();
      }
      else { // double band: plus/minus
        pm_band_integrand integrand_plus(bandstruct_,zeta_A);
        pm_band_integrand integrand_minus(bandstruct_,zeta_A, -1);
        I=twodsimpson(integrand_plus, bandstruct_.kx_min(), bandstruct_.ky_min(), bandstruct_.kx_max(), bandstruct_.ky_max(), L_);  //upper band
        I+=twodsimpson(integrand_minus, bandstruct_.kx_min(), bandstruct_.ky_min(), bandstruct_.kx_max(), bandstruct_.ky_max(), L_); // lower band
        I/=2.*bandstruct_.k_area();
      }
      G_omega_new(i,f)=I;
    }
  }
  //compute the new G0 
  for(unsigned f=0;f<G_omega.nflavor();++f)
    for(unsigned i=0;i<G_omega.nfreq();++i)
      G0_omega(i,f)=1./(Sigma(i,f)+1./G_omega_new(i,f));

  return G0_omega;
}


matsubara_green_function_t TwoDAFMHilbertTransformer::operator() (const matsubara_green_function_t & G_omega, 
                                                                 matsubara_green_function_t &G0_omega, 
                                                                 const double mu, const double h, const double beta)
{
  if(G_omega.nsite()!=1){throw std::logic_error("TwoDHilbertTransformer::operator(): don't know how to handle != 1 site.");}
  //if(G_omega.nflavor() !=2){throw std::logic_error("TwoDHilbertTransformer::operator(): don't know how to handle != 2 flavors.");}
  std::cout<<"AFM TwoD BS FS Hilbert Transformer"<<std::endl;
  matsubara_green_function_t Sigma(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
  matsubara_green_function_t G_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());

  for(unsigned i=0;i<G_omega.nfreq();++i){
    Sigma(i,0)=1./G0_omega(i,0)-1./G_omega(i,0);
    Sigma(i,1)=1./G0_omega(i,1)-1./G_omega(i,1);
  }
  //perform integration over bz
  for(unsigned i=0;i<G_omega.nfreq();++i){
    std::complex<double> iomegan(0,(2*i+1)*M_PI/beta);
    std::complex<double> zeta_A=iomegan + mu -h - Sigma(i,0);
    std::complex<double> zeta_B=iomegan + mu +h - Sigma(i,1);
    std::complex<double> I;
    afm_band_integrand integrand(bandstruct_,zeta_A,zeta_B);   // here the plus/minus band does give the same
    I=twodsimpson(integrand, bandstruct_.kx_min(), bandstruct_.ky_min(), bandstruct_.kx_max(), bandstruct_.ky_max(), L_) / bandstruct_.k_area();
    G_omega_new(i,0)=zeta_B*I;
    G_omega_new(i,1)=zeta_A*I;
  }
  //compute the new G0 
  for(unsigned i=0;i<G_omega.nfreq();++i){
    G0_omega(i,0)=1./(Sigma(i,0)+1./G_omega_new(i,0));
    G0_omega(i,1)=1./(Sigma(i,1)+1./G_omega_new(i,1));
  }
  return G0_omega;
}
