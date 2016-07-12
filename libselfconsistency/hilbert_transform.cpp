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
  p.define<double>("tprime", 0., "Next-nearest Hopping element for square lattice density of states ");
  p.define<double>("MU", "Chemical Potential");
  p.define<int>("L", 200, "Number of Brillouin zone points for square lattice integration ");
  p.define<std::string>("G0_OMEGA_INPUT", "Input filename for input G0");
  p.define<bool>("ATOMIC", false, "Generate atomic G0 to start from");
  p.define<double>("H_INIT", 0, "Initial (staggered) magnetic field");
  p.define<std::string>("sc.DOSFILE", "File containing table of Density of States");
  p.define<std::string>("sc.TWODBS", "Use Hilbert transform for two-dimensional band structure");

}
FrequencySpaceHilbertTransformer::FrequencySpaceHilbertTransformer(const alps::params &p):
            n_matsubara_(p["NMATSUBARA"]),
            n_orbital_(p["FLAVORS"]),
            beta_(p["BETA"]),
            mu_(p["MU"]),
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
  if(parms["sc.AFM_SELFCONSISTENCY"]){
    transform_ptr= new FSDOSHilbertTransformer(parms);
  }else{
    throw std::invalid_argument("put back the AFM DOS transformer, we have it in git.");
  }
  return transform_ptr;
}
FrequencySpaceHilbertTransformer *setup_2D_hilbert_transform(const alps::params &parms){
  FrequencySpaceHilbertTransformer *transform_ptr;
  if(parms["sc.SPIN_SYMMETRIZATION"] && parms["sc.AFM_SELFCONSISTENCY"]) throw std::invalid_argument("sc.AFM_SELFCONSISTENCY and sc.SPIN_SYMMETRIZATION? That doesn't make sense!");
  transform_ptr=new TwoDAFMHilbertTransformer(parms);
  return transform_ptr;
}
bool FrequencySpaceHilbertTransformer::read_initial_G0(single_site_matsubara_gf &G0_omega){
  {
    std::ifstream in_file(input_filename_.c_str()); if(!in_file.good()) return false;
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
void FrequencySpaceHilbertTransformer::set_c1_unity(single_site_matsubara_gf &G0_omega) const{
  typedef alps::gf::one_index_gf<double, alps::gf::index_mesh> density_matrix_type;
  density_matrix_type unity=density_matrix_type(alps::gf::index_mesh(2));
  unity.initialize();
  for(int i=0;i<G0_omega.mesh2().extent();++i)
    unity(alps::gf::index(i))=1;
  G0_omega.set_tail(1,unity);
}
void FrequencySpaceHilbertTransformer::set_c2(single_site_matsubara_gf &G0_omega) const{
  typedef alps::gf::one_index_gf<double, alps::gf::index_mesh> density_matrix_type;
  density_matrix_type c2=density_matrix_type(alps::gf::index_mesh(2));
  c2.initialize();
  for(int i=0;i<G0_omega.mesh2().extent();++i)
    c2(alps::gf::index(i))=(mu_-epsilon_av(i));
  G0_omega.set_tail(2,c2);
}

void FrequencySpaceHilbertTransformer::set_c3(single_site_matsubara_gf &G0_omega) const{
  typedef alps::gf::one_index_gf<double, alps::gf::index_mesh> density_matrix_type;
  density_matrix_type c3=density_matrix_type(alps::gf::index_mesh(2));
  c3.initialize();
  for(int i=0;i<G0_omega.mesh2().extent();++i)
    c3(alps::gf::index(i))=mu_*mu_-2.*mu_*epsilon_av(i)+epsilonsq_av(i);
  G0_omega.set_tail(3,c3);
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
        std::complex<double> zeta = iw+mu_;
        G0_omega(i, flavor) = 1./zeta;
      }
    }
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
    G0_omega=this->operator()(G_omega, G0_omega);
  }
  set_tails(G0_omega);
  write_initial_G0(G0_omega);
  return G0_omega;
}
single_site_matsubara_gf FSSemicircleHilbertTransformer::initial_G0()
{
  single_site_matsubara_gf G0_omega(alps::gf::omega_sigma_gf(alps::gf::matsubara_positive_mesh(beta_,n_matsubara_),alps::gf::index_mesh(n_orbital_)));
  if(read_initial_G0(G0_omega)==true){
    return G0_omega;
  }
  if (atomic_initial_G0_) {
    return FrequencySpaceHilbertTransformer::initial_G0();
  }
  std::cout<<"calculating initial G0_omega for semicircular DOS"<<std::endl;
  double tsq=t_*t_;
  for(alps::gf::matsubara_index i(0); i<G0_omega.mesh1().extent(); ++i) {
    std::complex<double> iw(0.,(2*i()+1)*M_PI/beta_);
    for(alps::gf::index flavor(0);flavor<G0_omega.mesh2().extent(); ++flavor) {
      std::complex<double> zeta = iw+mu_;
      std::complex<double> tmp = sqrt(zeta*zeta-4*tsq);
      tmp *= tmp.imag()<0 ? -1 : 1;
      G0_omega(i, flavor) = (zeta - tmp)/(2*tsq);
    }
  }
  set_tails(G0_omega);
  write_initial_G0(G0_omega);
  return G0_omega;
}
single_site_matsubara_gf FSSemicircleHilbertTransformer::operator()(const single_site_matsubara_gf& G_omega,single_site_matsubara_gf &G0_omega_ignored)
{
  std::cout<<"Single Site / Bethe Lattice Semi Circular Hilbert transform"<<std::endl;
  double tsq=t_*t_;
  single_site_matsubara_gf G0_omega(G0_omega_ignored);
  for(alps::gf::matsubara_index i(0); i<G0_omega.mesh1().extent(); ++i) {
    std::complex<double> iw(0.,(2*i()+1)*M_PI/beta_);
    for(alps::gf::index flavor(0);flavor<G0_omega.mesh2().extent(); flavor++) {
      G0_omega(i,flavor  ) =1./(iw + mu_ -tsq*G_omega(i,flavor));
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
    epsilon_.push_back(eps);
    dos_.push_back(d);
  }
  if(dos_.size()%2 !=0) throw std::invalid_argument("for Simpson precision: use a DOS with even number of integration points");
  normalize_dos();
  compute_moments();
}
void FSDOSHilbertTransformer::normalize_dos(){
  //normalize DOS to one:
  double s=0;
  for(unsigned int i=1;i<epsilon_.size()-2;i+=2){
    s+=4.*dos_[i]+2*dos_[i+1];
  }
  s+=dos_[0]+dos_[epsilon_.size()-1]+4.*dos_[epsilon_.size()-2];
  s/=(3.);
  std::cout<<"normalization constant: "<<s<<std::endl;
  std::cout<<"h is: "<<(epsilon_[1]-epsilon_[0])<<std::endl;
  for(unsigned int i=0;i<epsilon_.size()-1;++i){
    dos_[i]/=s;
  }
}
void FSDOSHilbertTransformer::compute_moments(){
  //find second moment of band structure:
  double epssqav=0;
  double epsav=0;
  for(unsigned int i=1;i<epsilon_.size()-2;i+=2){
    epssqav+=4.*(epsilon_[i]*epsilon_[i]*dos_[i])+2*(dos_[i+1]*epsilon_[i+1]*epsilon_[i+1]);
    epsav+=4.*(epsilon_[i]*dos_[i])+2*(dos_[i+1]*epsilon_[i+1]);
  }
  epssqav+=dos_[0]*epsilon_[0]*epsilon_[0]+dos_[epsilon_.size()-1]*epsilon_[epsilon_.size()-1]*epsilon_[epsilon_.size()-1]+4.*dos_[epsilon_.size()-2]*epsilon_[epsilon_.size()-2]*epsilon_[epsilon_.size()-2];
  epsav+=dos_[0]*epsilon_[0]+dos_[epsilon_.size()-1]*epsilon_[epsilon_.size()-1]+4.*dos_[epsilon_.size()-2]*epsilon_[epsilon_.size()-2];
  epssqav/=(3.);
  epsav/=(3.);

}


single_site_matsubara_gf FSDOSHilbertTransformer::operator()(const single_site_matsubara_gf &G_omega,single_site_matsubara_gf &G0_omega)
{
  //matsubara_green_function_t G0_omega(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor(), diagonal);
  single_site_matsubara_gf G_omega_new(G_omega);
  single_site_matsubara_gf G0_omega_new(G_omega);
  single_site_matsubara_gf sigma(G_omega);
  for(alps::gf::matsubara_index w(0); w<G0_omega.mesh1().extent(); ++w) {
    std::complex<double> iw(0.,(2*w()+1)*M_PI/beta_);
    for(alps::gf::index f(0);f<G0_omega.mesh2().extent(); ++f) {
      std::complex<double> g=0;
      sigma(w,f)= 1./G0_omega(w,f)-1./G_omega(w,f);
      std::complex<double> A=iw + mu_ -sigma(w, f);
      for(unsigned int n=0;n<dos_.size()-1;++n){
        g+=dos_[n]/(A-epsilon_[n]); //go back to higher order integrator!!
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

single_site_matsubara_gf TwoDAFMHilbertTransformer::operator()(const single_site_matsubara_gf & G_omega,single_site_matsubara_gf &G0_omega)
{
  //compute sigma
  single_site_matsubara_gf sigma(G_omega);
  single_site_matsubara_gf G(G_omega);
  if(paramagnet_){ spin_symmetrize(G0_omega); spin_symmetrize(G);}
  for(alps::gf::matsubara_index i(0); i<G0_omega.mesh1().extent(); ++i) {
    sigma(i,alps::gf::index(0))=1./G0_omega(i,alps::gf::index(0))-1./G(i,alps::gf::index(0));
    sigma(i,alps::gf::index(1))=1./G0_omega(i,alps::gf::index(1))-1./G(i,alps::gf::index(1));
  }
  //perform integration over bz
  for(alps::gf::matsubara_index i(0); i<G0_omega.mesh1().extent(); ++i) {
    G(i,alps::gf::index(0))=0.;
    G(i,alps::gf::index(1))=0.;
    std::complex<double> iomegan(0,(2*i()+1)*M_PI/beta_);
    std::complex<double> zeta_A=iomegan + mu_-sigma(i,alps::gf::index(0));
    std::complex<double> zeta_B=iomegan + mu_-sigma(i,alps::gf::index(1));
    afm_integrand integrand(t_,tprime_,zeta_A, zeta_B);

    std::complex<double> I=twodsimpson(integrand, -M_PI, -M_PI, M_PI, M_PI, L_);
    G(i,alps::gf::index(0))=1./(4*M_PI*M_PI)*zeta_B*I;
    G(i,alps::gf::index(1))=1./(4*M_PI*M_PI)*zeta_A*I;
  }
  //compute the new G0
  for(alps::gf::matsubara_index i(0); i<G0_omega.mesh1().extent(); ++i) {
    G0_omega(i,alps::gf::index(0))=1./(sigma(i,alps::gf::index(0))+1./G(i,alps::gf::index(0)));
    G0_omega(i,alps::gf::index(1))=1./(sigma(i,alps::gf::index(1))+1./G(i,alps::gf::index(1)));
  }
  if(paramagnet_){ spin_symmetrize(G0_omega);}
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

TwoDAFMHilbertTransformer::TwoDAFMHilbertTransformer(const alps::params& params):FrequencySpaceHilbertTransformer(params),
t_(params["t"]),
tprime_(params["tprime"]),
L_(params["L"]),
paramagnet_(params["sc.SPIN_SYMMETRIZATION"]){
  if(L_%2 !=0) throw std::runtime_error("please choose L to be even for Simpson integration for fast convergence!!");
}
