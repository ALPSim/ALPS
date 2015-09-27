/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2015 by Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/

#include<fstream>
#include<complex>
#include<cmath>
#include "hilbert_transform.h"

void FrequencySpaceHilbertTransformer::define_parameters(alps::params &p){
  p.define<double>("t", "Hopping element for semi circular & square lattice density of states ");
  p.define<double>("tprime", "Next-nearest Hopping element for square lattice density of states ");
  p.define<double>("MU", "Chemical Potential");
  p.define<int>("L", 200, "Number of Brillouin zone points for square lattice integration ");
  p.define<std::string>("G0_OMEGA_INPUT", "Input filename for input G0");
  p.define<bool>("ATOMIC", false, "Generate atomic G0 to start from");
  p.define<double>("H_INIT", 0, "Initial (staggered) magnetic field");
  p.define<std::string>("DOSFILE", "File contaning table of Density of States");
  p.define<std::string>("TWODBS", "Use Hilbert transform for two-dimensional band structure");

}
FrequencySpaceHilbertTransformer::FrequencySpaceHilbertTransformer(const alps::params &p):
        n_matsubara_(p["NMATSUBARA"]),
        n_orbital_(p["FLAVORS"]),
        beta_(p["BETA"]),
        mu_(p["MU"]),
        h_init_(p["H_INIT"]),
        input_filename_(p["G0_OMEGA_INPUT"].as<std::string>()),
        atomic_initial_G0_(p["ATOMIC"]){
}
FrequencySpaceHilbertTransformer *setup_Bethe_hilbert_transform(const alps::params &parms){
  std::cout<<"using Bethe lattice Hilbert transform"<<std::endl;
  return new FSSemicircleHilbertTransformer(parms);
}
FrequencySpaceHilbertTransformer *setup_DOS_hilbert_transform(const alps::params &parms){
  FrequencySpaceHilbertTransformer *transform_ptr;
  std::cout<<"using DOS Hilbert transform"<<std::endl;
  if(parms["AFM_SELFCONSISTENCY"]){
    transform_ptr= new FSDOSHilbertTransformer(parms);
  }else{
    transform_ptr= new AFM_FSDOSHilbertTransformer(parms);
  }
  return transform_ptr;
}
FrequencySpaceHilbertTransformer *setup_2D_hilbert_transform(const alps::params &parms){
  FrequencySpaceHilbertTransformer *transform_ptr;
  if(parms["ANTIFERROMAGNET"]==false){
    throw std::runtime_error("implement PM DOS integration!!");
  }else{
    if(parms["SPIN_SYMMETRIZATION"]==true){
      throw std::invalid_argument("ANTIFERROMAGNET and SPIN_SYMMETRIZATION? That doesn't make sense!");
    }
    transform_ptr=new TwoDAFMHilbertTransformer(parms);
  }
  return transform_ptr;
}
bool FrequencySpaceHilbertTransformer::read_initial_G0(single_site_matsubara_gf &G0_omega){
  {
    std::ifstream in_file(input_filename_); if(!in_file.good()) return false;
  }
  alps::hdf5::archive G0_input_file(input_filename_, alps::hdf5::archive::READ);
  std::cout<<"reading initial G0_omega"<<std::endl;
  G0_omega.load(G0_input_file, "/G0_omega");
  return true;
}
void FrequencySpaceHilbertTransformer::write_initial_G0(const single_site_matsubara_gf &G0_omega){
  alps::hdf5::archive G0_omega_file(input_filename_, alps::hdf5::archive::WRITE);
  G0_omega.save(G0_omega_file, "/G0_omega");
}
single_site_matsubara_gf FrequencySpaceHilbertTransformer::initial_G0(){
  single_site_matsubara_gf G0_omega(alps::gf::omega_sigma_gf(alps::gf::matsubara_positive_mesh(beta_,n_matsubara_),alps::gf::index_mesh(n_orbital_)));
  if(read_initial_G0(G0_omega)==true){
    return G0_omega;
  }
  if (atomic_initial_G0_) {
    std::cout<<"calculating initial G0_omega as atomic solution"<<std::endl;
    for(alps::gf::matsubara_index i(0); i<G0_omega.mesh1().extent(); ++i) {
      std::complex<double> iw(0.,(2*i()+1)*M_PI/beta_);
      for(alps::gf::index flavor(0);flavor<G0_omega.mesh2().extent(); ++flavor) {
        std::complex<double> zeta = iw+mu_+(flavor()%2==0 ? -h_init_ : h_init_);
        G0_omega(i, flavor) = 1./zeta;
      }
    }
    write_initial_G0(G0_omega);
    return G0_omega;
  }
  else {
    single_site_matsubara_gf G_omega(G0_omega);
    for(alps::gf::matsubara_index i(0); i<G0_omega.mesh1().extent(); ++i) {
      for(alps::gf::index flavor(0);flavor<G0_omega.mesh2().extent(); ++flavor) {
        //initializing G and G0 to 1 is a trick to get Sigma =0, which will produce a 'free' initial G0.
        G_omega(i,flavor)=1;
        G0_omega(i,flavor)=1;
      }
    }
    G0_omega=this->operator()(G_omega, G0_omega, mu_, h_init_, beta_);
    write_initial_G0(G0_omega);
  }
  return G0_omega;
}

single_site_matsubara_gf FSSemicircleHilbertTransformer::operator()(const single_site_matsubara_gf& G_omega,
    single_site_matsubara_gf &G0_omega_ignored,
    const double mu, const double h, const double beta)
{
  double tsq=t_*t_;
  single_site_matsubara_gf G0_omega(G_omega);
  for(alps::gf::matsubara_index i(0); i<G0_omega.mesh1().extent(); ++i) {
    std::complex<double> iw(0.,(2*i()+1)*M_PI/beta);
    for(alps::gf::index flavor(0);flavor<G0_omega.mesh2().extent(); flavor+=2) {
      G0_omega(i,flavor  ) =1./(iw + mu -h - tsq*G_omega(i,flavor+1));
      G0_omega(i,flavor+1) =1./(iw + mu +h - tsq*G_omega(i,flavor  ));
    }
  }
  return G0_omega;
}


FSDOSHilbertTransformer::FSDOSHilbertTransformer(const alps::params& params):FrequencySpaceHilbertTransformer(params)
{
  std::string dosfilename=params["DOSFILE"];
  std::ifstream dos_file(dosfilename.c_str());
  assert(dos_file.is_open());
  double eps, d;
  std::cout<<"using density of states: "<<std::endl;
  while(dos_file>>eps>>d){
    epsilon.push_back(eps);
    dos.push_back(d);
    std::cout<<eps<<" "<<d<<std::endl;
  }
  if(dos.size()%2 !=0) throw std::invalid_argument("for Simpson precision: use a DOS with even number of integration points");
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




single_site_matsubara_gf FSDOSHilbertTransformer::operator()(const single_site_matsubara_gf &G_omega,
    single_site_matsubara_gf &G0_omega,
    const double mu, const double h, const double beta)
{
  if(h!=0) throw std::invalid_argument(" we don't yet handle magnetic fields!");
  //matsubara_green_function_t G0_omega(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor(), diagonal);
  single_site_matsubara_gf G_omega_new(G_omega);
  single_site_matsubara_gf G0_omega_new(G_omega);
  single_site_matsubara_gf sigma(G_omega);
  for(alps::gf::matsubara_index w(0); w<G0_omega.mesh1().extent(); ++w) {
    std::complex<double> iw(0.,(2*w()+1)*M_PI/beta);
    for(alps::gf::index f(0);f<G0_omega.mesh2().extent(); ++f) {
      std::complex<double> g=0;
      sigma(w,f)= 1./G0_omega(w,f)-1./G_omega(w,f);
      std::complex<double> A=iw + mu -sigma(w, f);
      for(unsigned int n=0;n<dos.size()-1;++n){
        g+=dos[n]/(A-epsilon[n]); //go back to higher order integrator!!
      }
      G_omega_new(w, f)=g;
    }
  }
  for(alps::gf::matsubara_index w(0); w<G0_omega.mesh1().extent(); ++w) {
    for(alps::gf::index f(0);f<G0_omega.mesh2().extent(); ++f) {
      G0_omega_new(w,f)=1./(sigma(w,f)+1./G_omega_new(w,f));
    }
  }
  return G0_omega_new;
}




single_site_matsubara_gf AFM_FSDOSHilbertTransformer::operator()(const single_site_matsubara_gf &G_omega,
    single_site_matsubara_gf &G0_omega,
    const double mu, const double h, const double beta)
{
  single_site_matsubara_gf G_omega_new(G_omega);
  single_site_matsubara_gf G0_omega_new(G_omega);
  single_site_matsubara_gf sigma(G_omega);
  //compute self energy
  for(alps::gf::matsubara_index w(0); w<G0_omega.mesh1().extent(); ++w) {
    double wn=(2*w()+1.)*M_PI/beta;
    for(alps::gf::index f(0);f<G0_omega.mesh2().extent(); ++f) {
      sigma(w,f)=1./G0_omega(w,f)-1./G_omega(w,f);
    }
    std::complex<double>zeta_A=std::complex<double>(mu-h,wn)-sigma(w,alps::gf::index(0));
    std::complex<double>zeta_B=std::complex<double>(mu+h,wn)-sigma(w,alps::gf::index(1));

    //Simpson: integrate dos(e)/(zeta_A zeta_B-e^2)
    std::complex<double> I=0;
    for(unsigned int i=1;i<dos.size()-2;i+=2){
      I+=4.*integrand(i,dos,epsilon,zeta_A,zeta_B)+2.*integrand(i+1,dos,epsilon,zeta_A,zeta_B);
    }
    I+=integrand(0,dos,epsilon,zeta_A,zeta_B)+integrand(dos.size()-1,dos,epsilon,zeta_A,zeta_B)+4.*integrand(dos.size()-2,dos,epsilon,zeta_A,zeta_B);
    I/=3.;

    //compute the new G's:
    G_omega_new(w,alps::gf::index(0))=zeta_B*I; //formula 97 in review
    G_omega_new(w,alps::gf::index(1))=zeta_A*I;
    G0_omega_new(w,alps::gf::index(0))=1./(1./G_omega_new(w,alps::gf::index(0))+sigma(w,alps::gf::index(0)));
    G0_omega_new(w,alps::gf::index(1))=1./(1./G_omega_new(w,alps::gf::index(1))+sigma(w,alps::gf::index(1)));
  }

  return G0_omega_new;
}






single_site_matsubara_gf TwoDAFMHilbertTransformer::operator()(const single_site_matsubara_gf & G_omega,
    single_site_matsubara_gf &G0_omega,
    const double mu, const double h, const double beta)
{
  //compute sigma
  single_site_matsubara_gf sigma(G_omega);
  single_site_matsubara_gf G(G_omega);
  for(alps::gf::matsubara_index i(0); i<G0_omega.mesh1().extent(); ++i) {
    sigma(i,alps::gf::index(0))=1./G0_omega(i,alps::gf::index(0))-1./G_omega(i,alps::gf::index(0));
    sigma(i,alps::gf::index(1))=1./G0_omega(i,alps::gf::index(1))-1./G_omega(i,alps::gf::index(1));
  }
  //perform integration over bz
  for(alps::gf::matsubara_index i(0); i<G0_omega.mesh1().extent(); ++i) {
    G(i,alps::gf::index(0))=0.;
    G(i,alps::gf::index(1))=0.;
    std::complex<double> iomegan(0,(2*i()+1)*M_PI/beta);
    std::complex<double> zeta_A=iomegan + mu -h -sigma(i,alps::gf::index(0));
    std::complex<double> zeta_B=iomegan + mu +h -sigma(i,alps::gf::index(1));
    afm_integrand integrand(t_,tprime_,zeta_A, zeta_B);
    //std::cout<<i<<" "<<zeta_A<<" "<<zeta_B<<" "<<L_<<std::endl;
    std::complex<double> I=twodsimpson(integrand, -M_PI, -M_PI, M_PI, M_PI, L_);
    G(i,alps::gf::index(0))=1./(4*M_PI*M_PI)*zeta_B*I;
    G(i,alps::gf::index(1))=1./(4*M_PI*M_PI)*zeta_A*I;
  }
  //compute the new G0
  for(alps::gf::matsubara_index i(0); i<G0_omega.mesh1().extent(); ++i) {
    G0_omega(i,alps::gf::index(0))=1./(sigma(i,alps::gf::index(0))+1./G(i,alps::gf::index(0)));
    G0_omega(i,alps::gf::index(1))=1./(sigma(i,alps::gf::index(1))+1./G(i,alps::gf::index(1)));
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

