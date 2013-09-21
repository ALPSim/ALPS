/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2007 by Philipp Werner <werner@comp-phys.org>,
 *                            Matthias Troyer <troyer@comp-phys.org>,
 *                            Emanuel Gull <gullc@comp-phys.org>,
 *                            Sebastian Fuchs <fuchs@comp-phys.org>
 *               2012-2013 by Jakub Imriska <jimriska@phys.ethz.ch>
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


itime_green_function_t HilbertTransformer::initial_G0(const alps::Parameters& parms) const
{
  throw std::logic_error("not implemented - specify your Hilbert transformer");
}




itime_green_function_t SemicircleHilbertTransformer::operator()(const itime_green_function_t& G_tau, 
                                                                double mu, double h, double beta) const
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
  for(spin_t flavor=0;flavor<G_omega.nflavor(); flavor++){
    spin_t fbar=flavor%2==0?flavor+1:flavor-1;
    for(unsigned i=0; i<G_omega.nfreq(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      std::complex<double> zeta = iw + mu + (flavor%2 ? h : -h);
      G0_omega(i, flavor) = 1./(zeta - bethe_parms.second_moment(flavor)*G_omega(i,fbar));
    }
  }
  //std::cout<<"symmetrized G0 omega real: "<<std::endl;
  print_real_green_matsubara(std::cout, G0_omega, beta);
  fourier_ptr->backward_ft(G0_tau, G0_omega);
  return G0_tau;
}




itime_green_function_t SemicircleHilbertTransformer::initial_G0(const alps::Parameters& parms) const
{
  std::cout<<"SemicircleHilbertTransformer::initial_G0: ";
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
    GeneralFSHilbertTransformer hilbert(parms,false/*ignored*/);
    boost::shared_ptr<FourierTransformer> fourier_ptr;
    FourierTransformer::generate_transformer(parms, fourier_ptr);
    fourier_ptr->backward_ft(G0_tau, hilbert.initial_G0(parms));
  }
  
  if (parms.defined("G0TAU_input"))
    G0_tau.write((parms["G0TAU_input"]).c_str());
  
  return G0_tau;
}



GeneralFSHilbertTransformer::GeneralFSHilbertTransformer(const alps::Parameters& parms, bool ignored)
 : AFM(parms.value_or_default("ANTIFERROMAGNET",false)),
   bandstruct(BandstructureFactory(parms))
  {}

GeneralFSHilbertTransformer::GeneralFSHilbertTransformer(alps::Parameters& parms)
 : AFM(parms.value_or_default("ANTIFERROMAGNET",false)),
   bandstruct(BandstructureFactory(parms,true))
{
  bandstruct->set_parms(parms);
}


matsubara_green_function_t GeneralFSHilbertTransformer::initial_G0(const alps::Parameters& parms) const
{
  std::cout<<"GeneralFSHilbertTransformer::initial_G0: ";
  unsigned int n_matsubara=boost::lexical_cast<unsigned int>(parms["NMATSUBARA"]);
  unsigned int n_time=boost::lexical_cast<unsigned int>(parms["N"]);
  unsigned int n_orbital=parms.value_or_default("FLAVORS", 2);
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  double h = static_cast<double>(parms.value_or_default("H",0.));
  matsubara_green_function_t G0_omega(n_matsubara, n_orbital);

  if (parms.defined("INSULATING")) {
    std::cout<<"calculating insulating initial G0_omega"<<std::endl;
    for(unsigned int i=0; i<G0_omega.nfreq(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      for(spin_t flavor=0;flavor<G0_omega.nflavor(); flavor++) {
        std::complex<double> zeta = iw+mu+(flavor%2 ? h : -h);
        G0_omega(i, flavor) = 1./zeta;
      }
    }
  }
  else if (parms.defined("G0OMEGA_INPUT") && parms["G0OMEGA_INPUT"].length()>0) {
    std::cout<<"reading initial G0_omega"<<std::endl;
    std::ifstream check(parms["G0OMEGA_INPUT"].c_str());
    if(!check.is_open()) {
      std::cerr << "ERROR: could not open inital G0 file "<<parms["G0OMEGA_INPUT"]<<std::endl;
      throw std::runtime_error("GeneralFSHilbertTransformer::initial_G0: could not open inital G0 file");
    }
    else
      G0_omega.read(parms["G0OMEGA_INPUT"].c_str());
  }
  else {
    std::cout<<"calculating non-interacting initial G0_omega"<<std::endl;
    matsubara_green_function_t G_omega(n_matsubara, n_orbital);
    for(unsigned int i=0;i<n_matsubara;++i){
      for(unsigned int f=0;f<n_orbital;++f){
        G_omega(i,f)=1;
        G0_omega(i,f)=1;
      }
    }
    G0_omega=this->operator()(G_omega, G0_omega, mu, h, beta);
  }
  
  if (parms.defined("G0OMEGA_input"))   // it is not needed to store it by default, as it will be stored in the 1st iteration as G0_omega_1, G0_omegareal_1
    G0_omega.write((parms["G0OMEGA_input"]).c_str());

  return G0_omega;
}

matsubara_green_function_t GeneralFSHilbertTransformer::operator()(const matsubara_green_function_t& G_omega, 
                                                                      matsubara_green_function_t &G0_omega, 
                                                                      const double mu, const double h, const double beta) const
{
  if(G_omega.nsite()!=1){throw std::logic_error("GeneralFSHilbertTransformer::operator(): don't know how to handle systems with != 1 site.");}
  matsubara_green_function_t Sigma(1, 1, G_omega.nflavor());
  matsubara_green_function_t G_omega_new(1, 1, G_omega.nflavor());

  if (!AFM) {
    std::cout<<"GeneralFSHilbertTransformer: PM version; using: mu="<<mu<<", h="<<h<<", beta="<<beta<<std::endl;
    for(frequency_t w=0; w<G_omega.nfreq(); ++w){
      double wn=(2.*w+1)*M_PI/beta;
      for(spin_t f=0; f<G_omega.nflavor(); ++f){
        Sigma(0,f)=1./G0_omega(w,f)-1./G_omega(w,f);
        std::complex<double> zeta=std::complex<double>(mu-h,wn)-Sigma(0,f);
        
        G_omega_new(0,f)=bandstruct->HilbertIntegral_PM(zeta,f);
        G0_omega(w,f)=1./(1./G_omega_new(0,f)+Sigma(0,f));
      }
    }
  } else {
    std::cout<<"GeneralFSHilbertTransformer: AFM version; using: mu="<<mu<<", h="<<h<<", beta="<<beta<<std::endl;
    if(G_omega.nflavor()%2!=0){throw std::logic_error("GeneralFSHilbertTransformer::operator(): don't know how to handle odd number of flavors in the AFM case.");}
    for(spin_t f=0; f<G_omega.nflavor()/2; f+=2){
      for(frequency_t w=0; w<G_omega.nfreq(); ++w){
        Sigma(0,2*f)=1./G0_omega(w,2*f)-1./G_omega(w,2*f);
        Sigma(0,2*f+1)=1./G0_omega(w,2*f+1)-1./G_omega(w,2*f+1);
        double wn=(2.*w+1)*M_PI/beta;
        std::complex<double> zeta_0=std::complex<double>(mu-h,wn)-Sigma(0,2*f);
        std::complex<double> zeta_1=std::complex<double>(mu+h,wn)-Sigma(0,2*f+1);
        
        std::complex<double> integral=bandstruct->HilbertIntegral_AFM(zeta_0*zeta_1,2*f);
        
        //compute the new G's:
        G_omega_new(0,2*f)=zeta_1*integral; //formula 97 in review
        G_omega_new(0,2*f+1)=zeta_0*integral;
        G0_omega(w,2*f)=1./(1./G_omega_new(0,2*f)+Sigma(0,2*f));
        G0_omega(w,2*f+1)=1./(1./G_omega_new(0,2*f+1)+Sigma(0,2*f+1));
      }
    }
  }
  
  return G0_omega;
}



matsubara_green_function_t SemicircleFSHilbertTransformer::initial_G0(const alps::Parameters& parms) const {
  GeneralFSHilbertTransformer hilbert(parms,false/*ignored*/);
  /// NOTE: the initial_G0 calls the GeneralFSHilbertTransformer::operator(), which for Sigma=0 gives the G0
  /// NOTE: it WOULD NOT work with SemicircleFSHilbertTransformer::operator()
  return hilbert.initial_G0(parms);
}


matsubara_green_function_t SemicircleFSHilbertTransformer::operator()(const matsubara_green_function_t& G_omega, 
                                                                      matsubara_green_function_t &G0_omega_ignored, 
                                                                      const double mu, const double h, const double beta) const
{
  std::cout<<"SemicircleFSHilbertTransformer: using: mu:"<<mu<<" h: "<<h<<" beta: "<<beta<<" "<<G0_omega_ignored.nfreq()<<std::endl;
  matsubara_green_function_t G0_omega(G_omega);
  //formula according to review, p. 61, formula 221. 
  if(G_omega.nflavor()==1){ //special case. 
    for(unsigned i=0; i<G_omega.nfreq(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      G0_omega(i,0) =1./(iw + mu - bandstruct.second_moment(0)*G_omega(i,0));
    }
  }
  else{
    if(G_omega.nflavor()!=2){throw std::logic_error("SemicircleFSHilbertTransformer::operator(): don't know how to handle != 2 flavors in AFM case.");}
    for(unsigned flavor=0;flavor<G_omega.nflavor();flavor+=2){
      for(unsigned i=0; i<G_omega.nfreq(); i++) {
        std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
        G0_omega(i,flavor  ) =1./(iw + mu -h - bandstruct.second_moment(flavor)*G_omega(i,flavor+1)); 
        G0_omega(i,flavor+1) =1./(iw + mu +h - bandstruct.second_moment(flavor)*G_omega(i,flavor  )); 
      }
    }
  }
  return G0_omega;
}
