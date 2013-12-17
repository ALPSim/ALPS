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




itime_green_function_t HilbertTransformer::symmetrize(const itime_green_function_t& G_tau, const bool paramagnet) const
{
  itime_green_function_t G(G_tau);
  if (paramagnet) {
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
  if (parms.defined("ISOLATING")) {
    std::cout<<"calculating initial G0_tau"<<std::endl;
    int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
    int n_time=boost::lexical_cast<int>(parms["N"]);
    int n_flavor=parms.value_or_default("FLAVORS", 2);
    double beta = static_cast<double>(parms["BETA"]);
    double mu = static_cast<double>(parms["MU"]);
    double h = static_cast<double>(parms["H"]);
    matsubara_green_function_t G0_omega(n_matsubara,1, n_flavor, diagonal);
    itime_green_function_t G0_tau(n_time,1, n_flavor, diagonal);
    for(unsigned int i=0; i<G0_omega.nfreq(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      for(spin_t flavor=0;flavor<G0_omega.nflavor(); flavor++) {
        std::complex<double> zeta = iw+mu+(flavor%2 ? -h : h);
        G0_omega(i, flavor) = 1./zeta;
      }
    }
    G0_tau.write("G0_tau_input");
    G0_omega.write("G0_omega_input");
    return G0_tau;
  }
  else {
    std::cout<<"reading initial G0_tau"<<std::endl;
    int n_tau=boost::lexical_cast<int>(parms["N"]);
    int n_orbital=parms.value_or_default("FLAVORS", 2);
    itime_green_function_t G0_tau(n_tau+1, 1, n_orbital, diagonal);
    G0_tau.read(parms["G0TAU_INPUT"].c_str());
    return G0_tau;
  }
}




itime_green_function_t SemicircleHilbertTransformer::operator()(const itime_green_function_t& G_tau,
                                                                double mu, double h, double beta)
{
  //if(h!=0){ std::cout<<"ignoring magnetic field!"<<std::endl; abort(); }
  double tsq=t_*t_;
  matsubara_green_function_t G_omega(G_tau.ntime()-1, G_tau.nsite(), G_tau.nflavor(), diagonal);
  matsubara_green_function_t G0_omega(G_tau.ntime()-1, G_tau.nsite(), G_tau.nflavor(), diagonal);
  itime_green_function_t G0_tau(G_tau.ntime(), G_tau.nsite(), G_tau.nflavor(), diagonal);
  std::cerr << "SemicircleHilbertTransformer::operator(): Fouriertransformation of G_tau at this point is requested. Densities needed!\n";
  exit(1);
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  //FourierTransformer::generate_transformer(parms, fourier_ptr); ????
  fourier_ptr->time_to_frequency_ft(G_tau, G_omega);
  //std::cout<<"G omega real: "<<std::endl;
  print_real_green_matsubara(std::cout, G_omega, beta);
  for(unsigned int i=0; i<G_omega.nfreq(); i++) {
    std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
    for(spin_t flavor=0;flavor<G_omega.nflavor(); flavor++){
      std::complex<double> zeta = iw + mu + (flavor%2 ? -h : h);
      G0_omega(i, flavor) = 1./(zeta - tsq*G_omega(i, flavor%2==0?flavor+1:flavor-1));
    }
  }
  //std::cout<<"symmetrized G0 omega real: "<<std::endl;
  print_real_green_matsubara(std::cout, G0_omega, beta);
  fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega);
  return G0_tau;
}




itime_green_function_t SemicircleHilbertTransformer::initial_G0(const alps::Parameters& parms)
{
  if (parms.defined("G0TAU_INPUT") || parms.defined("ISOLATING"))
    return HilbertTransformer::initial_G0(parms);
  else {
    std::cout<<"calculating initial G0_tau"<<std::endl;
    int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
    int n_time=boost::lexical_cast<int>(parms["N"]);
    int n_flavor=parms.value_or_default("FLAVORS", 2);
    double beta = static_cast<double>(parms["BETA"]);
    double mu = static_cast<double>(parms["MU"]);
    double h = static_cast<double>(parms["H"]);
    //if(h!=0) {std::cout<<"we don't yet handle magnetic fields!"<<std::endl; abort(); }
    double tsq=t_*t_;
    matsubara_green_function_t G0_omega(n_matsubara,1,n_flavor, diagonal);
    itime_green_function_t G0_tau(n_time, 1, n_flavor, diagonal);
    for(unsigned int i=0; i<G0_omega.nfreq(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      for(spin_t flavor=0;flavor<G0_omega.nflavor(); flavor++) {
        std::complex<double> zeta = iw+mu+(flavor%2 ? -h : h);
        G0_omega(i, flavor) = (zeta - sqrt(zeta*zeta-4*tsq))/(2*tsq);
      }
    }
    G0_tau.write("G0_tau_input");
    G0_omega.write("G0_omega_input");
    return G0_tau;
  }
}




matsubara_green_function_t FrequencySpaceHilbertTransformer::initial_G0(const alps::Parameters& parms)
{
  //std::cout<<"General Hilbert: G0"<<std::endl;
  std::cout<<"reading initial G0_omega"<<std::endl;
  int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  int n_orbital=parms.value_or_default("FLAVORS", 2);
  matsubara_green_function_t G0_omega(n_matsubara, 1, n_orbital, diagonal);
  if (parms.defined("ISOLATING")) {
    std::cout<<"calculating initial G0_omega"<<std::endl;
    int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
    int n_time=boost::lexical_cast<int>(parms["N"]);
    int n_flavor=parms.value_or_default("FLAVORS", 2);
    double beta = static_cast<double>(parms["BETA"]);
    double mu = static_cast<double>(parms["MU"]);
    double h = static_cast<double>(parms["H"]);
    matsubara_green_function_t G0_omega(n_matsubara,1,n_flavor, diagonal);
    for(unsigned int i=0; i<G0_omega.nfreq(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      for(spin_t flavor=0;flavor<G0_omega.nflavor(); flavor++) {
        std::complex<double> zeta = iw+mu+(flavor%2 ? -h : h);
        G0_omega(i, flavor) = 1./zeta;
      }
    }
    G0_omega.write("G0_omega_input");
    boost::shared_ptr<FourierTransformer> fourier_ptr;
    FourierTransformer::generate_transformer(parms, fourier_ptr);
    itime_green_function_t G0_tau(n_time,1, n_flavor, diagonal);
    fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega);
    G0_tau.write("G0_tau_input");
    return G0_omega;
  }
  else {
    std::ifstream check(parms["G0OMEGA_INPUT"].c_str());
    if(!check.is_open()){
      std::cout<<"could not open inital G0 file "<<parms["G0OMEGA_INPUT"]<<" proceeding with free solution. "<<std::endl;
      matsubara_green_function_t G_omega(n_matsubara, 1, n_orbital, diagonal);
      for(int i=0;i<n_matsubara;++i){
        for(int f=0;f<n_orbital;++f){
          G_omega(i,f)=1;
          G0_omega(i,f)=1;
        }
      }
      double h     =parms.value_or_default("H",0.);
      double mu    =parms.value_or_default("MU",0.);
      int n_time=boost::lexical_cast<int>(parms["N"]);
      int n_flavor=parms.value_or_default("FLAVORS", 2);
      double beta  =parms.value_or_default("BETA",0.);
      G0_omega=this->operator()(G_omega, G0_omega, mu, h, beta);
      G0_omega.write("G0_omega_input");
      boost::shared_ptr<FourierTransformer> fourier_ptr;
      FourierTransformer::generate_transformer(parms, fourier_ptr);
      itime_green_function_t G0_tau(n_time,1,n_flavor, diagonal);
      fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega);
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
  std::cout<<"SCHT: using: mu:"<<mu<<" h: "<<h<<"beta: "<<beta<<" nomega: "<<G0_omega_ignored.nfreq()<<std::endl;
  double tsq=t_*t_;
  matsubara_green_function_t G0_omega(G_omega);
  if((G_omega.nflavor()%2)!=0)
    throw std::runtime_error("even number of flavors needed for AFM!");
  for(unsigned int flavor=0;flavor<G_omega.nflavor();flavor+=2){
    for(unsigned int i=0; i<G_omega.nfreq(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      G0_omega(i,flavor  ) =1./(iw + mu -h - tsq*G_omega(i,flavor+1));
      G0_omega(i,flavor+1) =1./(iw + mu +h - tsq*G_omega(i,flavor  ));
    }
  }
  return G0_omega;
}




matsubara_green_function_t FSSemicircleHilbertTransformer::initial_G0(const alps::Parameters& parms)
{
  //std::cout<<"Semicirc Hilbert: initial G0"<<std::endl;
  if ( parms.defined("ISOLATING"))
    return FrequencySpaceHilbertTransformer::initial_G0(parms);
  if (parms.defined("G0OMEGA_INPUT")){
    std::string name=parms["G0OMEGA_INPUT"];
    std::ifstream g0_file(name.c_str());
    if(g0_file.good())
      return FrequencySpaceHilbertTransformer::initial_G0(parms); //only open that file if it is good.
    else{
      std::cout<<"failed to open initial G0 file with name: "<<name<<std::endl;
    }
  }
  {
    std::cout<<"calculating initial G0_omega"<<std::endl;
    int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
    int n_time=boost::lexical_cast<int>(parms["N"]);
    int n_flavor=parms.value_or_default("FLAVORS", 2);
    double beta = static_cast<double>(parms["BETA"]);
    double mu = static_cast<double>(parms["MU"]);
    double h = static_cast<double>(parms["H"]);
    double tsq=t_*t_;
    matsubara_green_function_t G0_omega(n_matsubara,1,n_flavor, diagonal);
    for(unsigned int i=0; i<G0_omega.nfreq(); i++) {
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
    itime_green_function_t G0_tau(n_time,1,n_flavor, diagonal);
    fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega);
    G0_tau.write("G0_tau_input");
    //std::cout<<"written G0 in omega and tau, returning."<<std::endl;
    return G0_omega;
  }
}




FSDOSHilbertTransformer::FSDOSHilbertTransformer(const alps::Parameters& params)
{
  assert(params.defined("DOSFILE"));
  std::ifstream dos_file(params["DOSFILE"].c_str());
  assert(dos_file.is_open());
  double eps, d;
  std::cout<<"using density of states: "<<std::endl;
  while(dos_file>>eps>>d){
    epsilon.push_back(eps);
    dos.push_back(d);
    std::cout<<eps<<" "<<d<<std::endl;
  }
  if(dos.size()%2!=0){ std::cout<<"please use even number of DOS points"<<std::endl; abort(); }
  //normalize DOS to one:
  double s=0;
  for(unsigned int i=1;i<epsilon.size()-2;i+=2){
    s+=4.*dos[i]+2*dos[i+1];
  }
  s+=dos[0]+dos[epsilon.size()-1]+4.*dos[epsilon.size()-2];
  s/=(3.);
  std::cout<<"normalization constant: "<<s<<std::endl;
  std::cout<<"h is: "<<(epsilon[1]-epsilon[0])<<std::endl;
  for(unsigned int i=0;i<epsilon.size()-1;++i){
    dos[i]/=s;
  }
  double S=0;
  for(unsigned int i=0;i<epsilon.size()-1;++i){
    S+=dos[i];
  }
  std::cout<<"check: total sum is: "<<S<<std::endl;
  //find second moment of band structure:
  double epssqav=0;
  for(unsigned int i=1;i<epsilon.size()-2;i+=2){
    epssqav+=4.*(epsilon[i]*epsilon[i]*dos[i])+2*(dos[i+1]*epsilon[i+1]*epsilon[i+1]);
  }
  epssqav+=dos[0]*epsilon[0]*epsilon[0]+dos[epsilon.size()-1]*epsilon[epsilon.size()-1]*epsilon[epsilon.size()-1]+4.*dos[epsilon.size()-2]*epsilon[epsilon.size()-2]*epsilon[epsilon.size()-2];
  epssqav/=(3.);
  std::cout<<"second moment of band structure: "<<epssqav<<std::endl;
}




matsubara_green_function_t FSDOSHilbertTransformer::operator()(const matsubara_green_function_t &G_omega,
                                                               matsubara_green_function_t &G0_omega,
                                                               const double mu, const double h, const double beta)
{
  assert(G_omega.nsite()==1);
  if(h!=0){std::cout<<" we don't yet handle magnetic fields!"<<std::endl; abort(); }
  //matsubara_green_function_t G0_omega(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor(), diagonal);
  matsubara_green_function_t G_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor(), diagonal);
  matsubara_green_function_t G0_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor(), diagonal);
  matsubara_green_function_t sigma(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor(), diagonal);
  for(spin_t f=0;f<G_omega.nflavor();++f){
    //std::cout<<"check of G0omega: "<<f<<" "<<G0_omega(0,f)<<" "<<G_omega(0,f);
    for(frequency_t w=0;w<G_omega.nfreq();++w){
      std::complex<double> g=0;
      sigma(w,f)= 1./G0_omega(w,f)-1./G_omega(w,f);
      std::complex<double> A=std::complex<double>(0, (2*w+1)*M_PI/beta) + mu -sigma(w, f);
      for(unsigned int n=0;n<dos.size()-1;++n){
        //double delta_epsilon = epsilon[n]-epsilon[n+1];
        //double delta_dos     = dos[n]    -dos    [n+1];
        //perform analytic integration of D(e)/(A-e) for D known at points
        //e0, e1 at the end of the interval.
        //g+=1./(delta_epsilon)*(delta_dos*delta_epsilon+(A*delta_dos +dos[n+1]*epsilon[n]-dos[n]*epsilon[n+1])
        //                       *(log(A-epsilon[n]) - log(A-epsilon[n+1])));
        //g+=0.5*(dos[n+1]-dos[n])*(log(A-epsilon[n]) + log(A-epsilon[n+1]));
        g+=dos[n]/(A-epsilon[n]); //go back to higher order integrator!!
      }
      G_omega_new(w, f)=g;
    }
    //std::cout<<"G omega new is: "<<G_omega_new(0,f)<<std::endl;
  }
  if(G_omega.nflavor() !=2){std::cout<<"don't know how to handle != 2 flavors."<<std::endl; abort();}
  for(spin_t f=0;f<G_omega.nflavor();++f){
    for(frequency_t w=0;w<G_omega.nfreq();++w){
      G0_omega_new(w,f)=1./(sigma(w,f)+1./G_omega_new(w,f));
    }
  }
  //std::cout<<"G0_omega_new"<<" 0 "<<G0_omega_new(0,0)<<std::endl;
  //std::cout<<"G0_omega_new"<<" 1 "<<G0_omega_new(0,1)<<std::endl;
  return G0_omega_new;
}




matsubara_green_function_t AFM_FSDOSHilbertTransformer::operator()(const matsubara_green_function_t &G_omega,
                                                                   matsubara_green_function_t &G0_omega,
                                                                   const double mu, const double h, const double beta)
{
  if(G_omega.nsite()!=1){ std::cout<<"don't know how to handle n_site >1. Use Sebastian's cluster loop for that. aborting."<<std::endl; abort();}
  if(G_omega.nflavor()!=2){ std::cout<<"don't know how to handle n_flavor!=2. aborting."<<std::endl; abort();}
  matsubara_green_function_t Sigma(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor(), diagonal);
  matsubara_green_function_t G_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor(), diagonal);
  matsubara_green_function_t G0_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor(), diagonal);
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
    if(dos.size()%2 !=0){std::cout<<"for Simpson precision: use a DOS with even number of integration points. "<<std::endl; abort(); }
    std::complex<double> I=0;
    for(unsigned int i=1;i<dos.size()-2;i+=2){
      I+=4.*integrand(i,dos,epsilon,zeta_A,zeta_B)+2.*integrand(i+1,dos,epsilon,zeta_A,zeta_B);
    }
    I+=integrand(0,dos,epsilon,zeta_A,zeta_B)+integrand(dos.size()-1,dos,epsilon,zeta_A,zeta_B)+4.*integrand(dos.size()-2,dos,epsilon,zeta_A,zeta_B);
    I/=3.;
    /*for(int i=0;i<dos.size();++i){
     I+=dos[i]/(zeta_A*zeta_B - epsilon[i]*epsilon[i]);
     }*/
    //compute the new G's:
    G_omega_new(w,0)=zeta_B*I; //formula 97 in review
    G_omega_new(w,1)=zeta_A*I;
    G0_omega_new(w,0)=1./(1./G_omega_new(w,0)+Sigma(w,0));
    G0_omega_new(w,1)=1./(1./G_omega_new(w,1)+Sigma(w,1));
    //std::cout<<"G0: "<<G0_omega_new(w,0)<<" "<<G0_omega_new(w,1)<<std::endl;
  }
  //std::cout<<std::endl<<std::endl;
  /*std::cout<<"sigma: "<<std::endl;
   for(int i=0;i<Sigma.nfreq();++i){
   std::cout<<(2.*i+1)*M_PI/beta<<" "<<Sigma(i,0).real()<<" "<<Sigma(i,0).imag()<<" "<<Sigma(i,1).real()<<" "<<Sigma(i,1).imag()<<std::endl;
   }
   std::cout<<"g0 omega new: "<<std::endl;
   for(int i=0;i<G0_omega_new.nfreq();++i){
   std::cout<<(2.*i+1)*M_PI/beta<<" "<<G0_omega_new(i,0).real()<<" "<<G0_omega_new(i,0).imag()<<" "<<G0_omega_new(i,1).real()<<" "<<G0_omega_new(i,1).imag()<<std::endl;
   }*/
  return G0_omega_new;
}






matsubara_green_function_t TwoDAFMHilbertTransformer::operator()(const matsubara_green_function_t & G_omega,
                                                                 matsubara_green_function_t &G0_omega,
                                                                 const double mu, const double h, const double beta)
{
  //compute sigma
  matsubara_green_function_t sigma(G_omega);
  matsubara_green_function_t G(G_omega);
  for(unsigned int i=0;i<G_omega.nfreq();++i){
    sigma(i,0)=1./G0_omega(i,0)-1./G_omega(i,0);
    sigma(i,1)=1./G0_omega(i,1)-1./G_omega(i,1);
  }
  //perform integration over bz
  for(unsigned int i=0;i<G_omega.nfreq();++i){
    G(i,0)=0.;
    G(i,1)=0.;
    std::complex<double> iomegan(0,(2*i+1)*M_PI/beta);
    std::complex<double> zeta_A=iomegan + mu -h -sigma(i,0);
    std::complex<double> zeta_B=iomegan + mu +h -sigma(i,1);
    afm_integrand integrand(t_,tprime_,zeta_A, zeta_B);
    //std::cout<<i<<" "<<zeta_A<<" "<<zeta_B<<" "<<L_<<std::endl;
    std::complex<double> I=twodsimpson(integrand, -M_PI, -M_PI, M_PI, M_PI, L_);
    G(i,0)=1./(4*M_PI*M_PI)*zeta_B*I;
    G(i,1)=1./(4*M_PI*M_PI)*zeta_A*I;
  }
  //compute the new G0
  for(unsigned int i=0;i<G_omega.nfreq();++i){
    G0_omega(i,0)=1./(sigma(i,0)+1./G(i,0));
    G0_omega(i,1)=1./(sigma(i,1)+1./G(i,1));
  }
  return G0_omega;
}




#ifdef IGOR

extern "C" void kintegration_(const int *nr, const double *W, const double *mu,
                              const int *n_band, const void *sig_qmc, void *Gf_qmc);



matsubara_green_function_t IgorsHilbertTransformer::operator()(const matsubara_green_function_t& G_omega,
                                                               matsubara_green_function_t &G0_omega,
                                                               const double /*mu*/, const double /*h*/, const double beta)
{
  assert(G_omega.nflavor()%2==0);
  assert(G_omega.nsite()==1);
  int nfreq = G_omega.nfreq();
  int nband = G_omega.nflavor()/2;
  std::vector<std::complex<double> > selfenergy_up(nfreq*nband);
  std::vector<std::complex<double> > selfenergy_down(nfreq*nband);
  for(unsigned int flavor=0; flavor<G_omega.nflavor(); flavor+=2) {
    for(unsigned int i=0; i<nfreq; i++) {
      selfenergy_up[i*nband+flavor/2] = 1./G0_omega(i, flavor) - 1./G_omega(i, flavor);
      selfenergy_down[i*nband+flavor/2] = 1./G0_omega(i, flavor+1) - 1./G_omega(i, flavor+1);
    }
  }
  double omega_n[nfreq];
  for(unsigned int n=0; n<nfreq; n++)
    omega_n[n] = (2*n+1)*M_PI/beta;
    std::vector<std::complex<double> > G0_omega_new_up(nfreq*nband);
    std::vector<std::complex<double> > G0_omega_new_down(nfreq*nband);
    kintegration_(&nfreq, omega_n, &mu, &nband, &selfenergy_up[0], &G0_omega_new_up[0]);
    kintegration_(&nfreq, omega_n, &mu, &nband, &selfenergy_down[0], &G0_omega_new_down[0]);
    for(unsigned int flavor=0; flavor<nband; ++flavor) {
      for(unsigned int i=0; i<nfreq; i++) {
        G0_omega_new_up[i*nband+flavor] = 1./(1./G0_omega_new_up[i*nband+flavor] + selfenergy_up[i*nband+flavor]);
        G0_omega_new_down[i*nband+flavor] = 1./(1./G0_omega_new_down[i*nband+flavor] + selfenergy_down[i*nband+flavor]);
      }
    }
  matsubara_green_function_t G0_omega_new(G_omega);
  for(unsigned int flavor=0; flavor<G0_omega_new.nflavor(); flavor+=2) {
    for(unsigned int i=0; i<nfreq; i++) {
      G0_omega_new(i, flavor) = G0_omega_new_up[i*nband+flavor/2];
      G0_omega_new(i, flavor+1) = G0_omega_new_down[i*nband+flavor/2];
    }
  }
  return G0_omega_new;
}




matsubara_green_function_t IgorsHilbertTransformer::initial_G0(const alps::Parameters& parms)
{
  if (parms.defined("G0OMEGA_INPUT"))
    return FrequencySpaceHilbertTransformer::initial_G0(parms);
  else {
    std::cout<<"calculating initial G0_omega"<<std::endl;
    int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
    //int n_time=boost::lexical_cast<int>(parms["N"]);
    int n_flavor=parms.value_or_default("FLAVORS", 2);
    double beta = static_cast<double>(parms["BETA"]);
    double mu = static_cast<double>(parms["MU"]);
    double h = static_cast<double>(parms["H"]);
    matsubara_green_function_t G0_omega(n_matsubara,n_flavor);
    for(unsigned int i=0; i<G0_omega.nfreq(); i++) {
      for(spin_t flavor=0;flavor<G0_omega.nflavor(); flavor++) {
        G0_omega(i, flavor) = 1;
      }
    }
    matsubara_green_function_t G_omega(G0_omega);
    matsubara_green_function_t G0_omega_new = operator()(G_omega, G0_omega, mu, h, beta);
    return G0_omega_new;
  }
}

#endif
