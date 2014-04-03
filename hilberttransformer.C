/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2014 by Emanuel Gull <gullc@comp-phys.org>
 *
 * To be used with permission by the authors
 *
 *****************************************************************************/

/* $Id: hilberttransformer.C 338 2009-01-20 01:22:05Z fuchs $ */



#include "hilberttransformer.h"
#include <functional>


FrequencySpaceHilbertTransformer *setup_Bethe_hilbert_transform(const alps::Parameters &parms){
  std::cout<<"using Bethe lattice Hilbert transform"<<std::endl;
  return new FSSemicircleHilbertTransformer(boost::lexical_cast<double>(parms["t"]));
}
FrequencySpaceHilbertTransformer *setup_DOS_hilbert_transform(const alps::Parameters &parms){
  FrequencySpaceHilbertTransformer *transform_ptr;
  std::cout<<"using DOS Hilbert transform"<<std::endl;
  if(parms.value_or_default("ANTIFERROMAGNET",false)){
    transform_ptr= new FSDOSHilbertTransformer(parms);
  }else{
    transform_ptr= new AFM_FSDOSHilbertTransformer(parms);
  }
  return transform_ptr;
}
FrequencySpaceHilbertTransformer *setup_2D_hilbert_transform(const alps::Parameters &parms){
  FrequencySpaceHilbertTransformer *transform_ptr;
  if(!parms.defined("ANTIFERROMAGNET") || ((bool)(parms.value_or_default("ANTIFERROMAGNET",false))==false)){
    if((bool)(parms.value_or_default("PARAMAGNET",true))==false){
      std::cerr<<"AFM and PM? redundant parameters! set paramagnet to false."<<std::endl; abort();
    }
    std::cerr<<"implement PM DOS integration!!"<<std::endl;
    abort();
  }else{
    if((bool)(parms["PARAMAGNET"])==true){
      std::cerr<<"AFM and PM? redundant parameters! set ANTIFERROMAGNET to false."<<std::endl; abort();
    }
    transform_ptr=new TwoDAFMHilbertTransformer(parms);
  }
  return transform_ptr;
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
    double h = static_cast<double>(parms["H_INIT"]);
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
      double h     =parms.value_or_default("H_INIT",0.);
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


std::complex<double> TwoDAFMHilbertTransformer::twodsimpson(const integrand &f, double ax, double ay, double bx, double by, int N){
  std::complex<double>       result=0.;
  double       h=(bx-ax)/(2.*N);
  double       k=(by-ay)/(2.*N);

  result =  f(ax,ay) + f(ax,by) + f(bx, ay) +f(bx,by) ;

  // values between boundaries
  for ( int i = 1; i <= N; ++i ) {
    result += 4.*f(ax,ay+k*(2*i-1));
  }
  for ( int i = 1; i <= N-1; ++i ) {
    result += 2.*f(ax,ay+k*2*i);
  }
  for ( int i = 1; i <= N; ++i ) {
    result += 4.*f(bx,ay+k*(2*i-1));
  }
  for ( int i = 1; i <= N-1; ++i ) {
    result += 2.*f(bx,ay+k*2*i);
  }


  for ( int i = 1; i <= N; ++i ) {
    result += 4.*f(ax+h*(2*i-1),ay);
  }
  for ( int i = 1; i <= N-1; ++i ) {
    result += 2.*f(ax+h*(2*i),ay);
  }
  for ( int i = 1; i <= N; ++i ) {
    result += 4.*f(ax+h*(2*i-1),by);
  }
  for ( int i = 1; i <= N-1; ++i ) {
    result += 2.*f(ax+h*(2*i),by);
  }

  //inner part
  //
  for(int i=1;i<=N;++i){
    for(int j=1;j<=N;++j){
      result+=16.*f(ax+h*(2*i-1), ay+k*(2*j-1));
    }
  }
  for(int i=1;i<N;++i){
    for(int j=1;j<=N;++j){
      result+=8.*f(ax+h*(2*i-1), ay+k*(2*j));
    }
  }
  for(int i=1;i<=N;++i){
    for(int j=1;j<N;++j){
      result+=8.*f(ax+h*(2*i), ay+k*(2*j-1));
    }
  }
  for(int i=1;i<N;++i){
    for(int j=1;j<N;++j){
      result+=4.*f(ax+h*(2*i), ay+k*(2*j));
    }
  }

  result *= h*k/9.;
  return result;
}

