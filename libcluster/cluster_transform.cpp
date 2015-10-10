/*****************************************************************************
 *
 * Copyright (C) 2007-2015 by Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/

#include <iostream>
#include <string>

#include "alps/numeric/vector_functions.hpp"
#include "alps/hdf5/vector.hpp"
#include "alps/gf/fourier.hpp"
#include "alps/params.hpp"

#include "cluster_transform.h"

void define_dca_parameters(alps::params &p){
  ClusterHelper::define_parameters(p);
  p.define<double>("dca.t_flavor", 0, "Hopping for mass imbalance for second orbital");
  p.define<double>("dca.t", "nearest-neighbor hopping parameter and general unit");
  p.define<double>("dca.tprime", 0, "next-nearest neighbor hopping parameter");
  p.define<double>("dca.H_INIT", 0, "Initial (staggered) magnetic field");
  p.define<double>("dca.H", 0, "magnetic field");
  p.define<double>("dca.MU", "Chemical Potential");
  p.define<double>("dca.BETA", "Inverse temperature");
  p.define<bool>("dca.IMBALANCE", false, "set to true for mass imbalance simulation");
  p.define<bool>("dca.PHASE", false, "set to true for simulation with phase");
  p.define<double>("dca.PHASE_0", 0, "Phase in direction 0 ");
  p.define<double>("dca.PHASE_1", 0, "Phase in direction 1 ");
  p.define<double>("dca.PHASE_2", 0, "Phase in direction 2 ");
  p.define<int>("dca.NMATSUBARA", "Number fo Matsubara frequencies");
  p.define<int>("dca.N", "Number of imaginary time points");
  p.define<int>("dca.FLAVORS", "Number of orbitals (usually 2)");
  p.define<std::string>("dca.G0_OMEGA_INPUT", "Location of the input G0 file");
  p.define<std::string>("dca.INITIAL_G0", "FREE", "choose ISOLATED or FREE for initial G0s");
}
NormalStateClusterTransformer::NormalStateClusterTransformer(const alps::params& params) : ClusterHelper(params),
    mu_(params["dca.MU"]),
    h_(params["dca.H"]),
    half_filling_(0.),
    t_(params["dca.t"]),
    t2_(params["dca.tprime"]),
    t_flavor_(params["dca.t_flavor"]),
    n_flavors_(params["dca.FLAVORS"]),
    use_phase_(params["dca.PHASE"]),
    use_imbalance_(params["dca.IMBALANCE"]),
    epsilon_kckl_(boost::extents[n_cluster_sites()][n_lattice_momenta()][n_flavors_]),
    weight_kl_(n_lattice_momenta()),
    nfreq_(params["dca.NMATSUBARA"]),
    beta_(params["dca.BETA"]),
    g0_input_path_(params["dca.G0_OMEGA_INPUT"].as<std::string>()),
    initial_g0_(params["dca.INITIAL_G0"].as<std::string>())
{
  compute_integration_weight_and_dispersion();
  compute_dispersion_weights();
  if(use_phase_) setup_phase_vector(params);
}

bool NormalStateClusterTransformer::read_initial_G0(cluster_matsubara_kspace_gf &G0_omega) const{
  std::cout<<"reading initial G0_omega from: "<<g0_input_path_<<std::endl;
  {
    std::ifstream in_file(g0_input_path_);
    if(!in_file.good())
      std::cout<<"WARNING: Input file could not be found, computing initial solution."<<std::endl;
    return false;
  }

  alps::hdf5::archive ar(g0_input_path_);
  G0_omega.load(ar, "/G0");
  return true;
}
void NormalStateClusterTransformer::compute_isolated_cluster_G0(cluster_matsubara_kspace_gf &G0_omega) const{
  std::cout << "calculating G0_omega for isolated cluster"<<std::endl;
  std::vector<double> k(dimension());
  for (alps::gf::matsubara_index freq(0); freq<nfreq_; ++freq) {
    const std::complex<double> iomega_n(0, G0_omega.mesh1().points()[freq()]);
    for (alps::gf::momentum_index kc(0); kc<n_cluster_sites(); ++kc) {
      for (unsigned int d=0; d<dimension(); ++d)
        k[d] = cluster_momentum(kc(),d);
      const std::complex<double> G0 = use_phase_?1./(iomega_n+mu_-epsilon_phi(k,phi_)):1./(iomega_n+mu_-epsilon(k));
      for (alps::gf::index flavor(0); flavor<n_flavors_; ++flavor) {
        G0_omega(freq, kc, flavor) = G0;
      }
    }
  }
}
void NormalStateClusterTransformer::compute_noninteracting_G0(cluster_matsubara_kspace_gf &G0_omega) const{
  cluster_matsubara_kspace_gf zero_selfenergy(G0_omega);
  zero_selfenergy.initialize();
  matsubara_k_sum(G0_omega, zero_selfenergy);
}
cluster_matsubara_kspace_gf NormalStateClusterTransformer::initial_G0()const{

  cluster_matsubara_kspace_gf G0_omega(alps::gf::omega_k_sigma_gf(
      alps::gf::matsubara_positive_mesh(beta_,nfreq_),
      momentum_mesh(),
      alps::gf::index_mesh(n_flavors_))
  );

  if(read_initial_G0(G0_omega)) return G0_omega;

  if(initial_g0_=="FREE")
    compute_noninteracting_G0(G0_omega);
  else if (initial_g0_=="ISOLATED")
    compute_isolated_cluster_G0(G0_omega);
  else
    throw std::invalid_argument("choose FREE or ISOLATED as initial G0 parameters or provide an input hdf5 file");

  {
    alps::hdf5::archive ar(g0_input_path_, alps::hdf5::archive::WRITE);
    G0_omega.save(ar, "/G0");

    std::ofstream G0_omega_input_file("G0_omega_input.dat"); G0_omega_input_file<<G0_omega;
  }
  {
    cluster_itime_kspace_gf G0_tau(alps::gf::itime_k_sigma_gf(
        alps::gf::itime_mesh(beta_,nfreq_),
        momentum_mesh(),
        alps::gf::index_mesh(n_flavors_))
    );
    alps::gf::fourier_frequency_to_time(G0_omega, G0_tau);
    std::ofstream G0_tau_input_file("G0_tau_input.dat"); G0_tau_input_file<<G0_tau;
  }
  return G0_omega;
}

alps::gf::momentum_index_mesh NormalStateClusterTransformer::momentum_mesh() const{
  alps::gf::momentum_index_mesh::container_type points(boost::extents[n_cluster_sites()][dimension()]);
  for(int i=0;i<n_cluster_sites();++i){
    for(int j=0;j<dimension();++j){
      points[i][j]=cluster_momentum(i,j);
    }
  }
  return alps::gf::momentum_index_mesh(points);
}
alps::gf::real_space_index_mesh NormalStateClusterTransformer::real_space_mesh() const{
  alps::gf::real_space_index_mesh::container_type points(boost::extents[n_cluster_sites()][dimension()]);
  for(int i=0;i<n_cluster_sites();++i){
    for(int j=0;j<dimension();++j){
      points[i][j]=cluster_coordinate(i,j);
    }
  }
  return alps::gf::real_space_index_mesh(points);
}

void NormalStateClusterTransformer::spin_symmetrize(cluster_itime_kspace_gf& gf) const{
  if(gf.mesh3().extent()!=2) throw std::runtime_error("mesh should only have two flavors.");
  for (alps::gf::itime_index om=alps::gf::itime_index(0); om<gf.mesh1().extent(); ++om) {
    for (alps::gf::momentum_index ii=alps::gf::momentum_index(0); ii<gf.mesh2().extent(); ++ii) {
      double mean=(gf(om,ii,alps::gf::index(0))+gf(om,ii,alps::gf::index(1)))/2.;
      gf(om,ii,alps::gf::index(0))=mean;
      gf(om,ii,alps::gf::index(1))=mean;
    }
  }
}
///perform spin symmetrization of matsubara k-space Green's function
void NormalStateClusterTransformer::spin_symmetrize(cluster_matsubara_kspace_gf& gf) const{
  if(gf.mesh3().extent()!=2) throw std::runtime_error("mesh should only have two flavors.");
  for (alps::gf::matsubara_index om=alps::gf::matsubara_index(0); om<gf.mesh1().extent(); ++om) {
    for (alps::gf::momentum_index ii=alps::gf::momentum_index(0); ii<gf.mesh2().extent(); ++ii) {
      std::complex<double> mean=(gf(om,ii,alps::gf::index(0))+gf(om,ii,alps::gf::index(1)))/2.;
      gf(om,ii,alps::gf::index(0))=mean;
      gf(om,ii,alps::gf::index(1))=mean;
    }
  }
}
///perform spin symmetrization of imaginary time k-space cluster Green's function
void NormalStateClusterTransformer::spin_symmetrize(cluster_itime_realspace_gf& gf) const{
  if(gf.mesh3().extent()!=2) throw std::runtime_error("mesh should only have two flavors.");
  for (alps::gf::itime_index om=alps::gf::itime_index(0); om<gf.mesh1().extent(); ++om) {
    for (alps::gf::real_space_index ii=alps::gf::real_space_index(0); ii<gf.mesh2().extent(); ++ii) {
      for (alps::gf::real_space_index jj=alps::gf::real_space_index(0); jj<gf.mesh3().extent(); ++jj) {
        double mean=(gf(om,ii,jj,alps::gf::index(0))+gf(om,ii,jj,alps::gf::index(1)))/2.;
        gf(om,ii,jj,alps::gf::index(0))=mean;
        gf(om,ii,jj,alps::gf::index(1))=mean;
      }
    }
  }
}
///perform spin symmetrization of matsubara k-space Green's function
void NormalStateClusterTransformer::spin_symmetrize(cluster_matsubara_realspace_gf& gf) const{
  if(gf.mesh3().extent()!=2) throw std::runtime_error("mesh should only have two flavors.");
  for (alps::gf::matsubara_index om=alps::gf::matsubara_index(0); om<gf.mesh1().extent(); ++om) {
    for (alps::gf::real_space_index ii=alps::gf::real_space_index(0); ii<gf.mesh2().extent(); ++ii) {
      for (alps::gf::real_space_index jj=alps::gf::real_space_index(0); jj<gf.mesh3().extent(); ++jj) {
        std::complex<double> mean=(gf(om,ii,jj,alps::gf::index(0))+gf(om,ii,jj,alps::gf::index(1)))/2.;
        gf(om,ii,jj,alps::gf::index(0))=mean;
        gf(om,ii,jj,alps::gf::index(1))=mean;
      }
    }
  }
}
struct Forward {
  template<class GF_RS, class GF_KS> static void add(GF_KS & G_new, const GF_RS& G,
      const int i, const int j, const int k,
      double factor) {
    for (typename GF_KS::mesh1_type::index_type t(0); t<G.mesh1().extent(); ++t) {
      for (alps::gf::index flavor(0); flavor<G.mesh4().extent(); ++flavor) {
        G_new(t, alps::gf::momentum_index(k), flavor) +=
            factor * G(t, alps::gf::real_space_index(i), alps::gf::real_space_index(j), flavor);
      }
    }
  }
  template<class GF_RS, class GF_KS>
  static void add(GF_RS& G_new, const GF_KS& G,
      const int t, const int i, const int j, const int k,
      double factor) {
    for(alps::gf::index s1(0); s1<2;++s1){
      for(alps::gf::index s2(0); s2<2;++s2){
        G_new(GF_KS::mesh1_type::index_type(t), alps::gf::momentum_index(k), s1,s2) +=
            factor * G(GF_RS::mesh1_type::index_type(t), alps::gf::real_space_index(i), alps::gf::real_space_index(j), s1, s2);
      }
    }
  }
};



struct Backward {
  template<class GF_RS, class GF_KS> static void add(GF_RS & G_new, const GF_KS& G,
      const int i, const int j, const int k,
      double factor) {
    for (typename GF_RS::mesh1_type::index_type t(0); t<G.mesh1().extent(); ++t) {
      for (alps::gf::index flavor(0); flavor<G.mesh3().extent(); ++flavor) {
        G_new(t, alps::gf::real_space_index(i), alps::gf::real_space_index(j), flavor) +=
            factor * G(t, alps::gf::momentum_index(k), flavor);
      }
    }
  }
  template<class GF_RS, class GF_KS>
  static void add(GF_RS& G_new, const GF_KS& G,
      const int t, const int i, const int j, const int k,
      double factor) {
    for(alps::gf::index s1(0); s1<2;++s1){
      for(alps::gf::index s2(0); s2<2;++s2){
        G_new(GF_RS::mesh1_type::index_type(t), alps::gf::real_space_index(i), alps::gf::real_space_index(j), s1,s2) +=
            factor * G(GF_RS::mesh1_type::index_type(t), alps::gf::momentum_index(k), s1, s2);
      }
    }
  }
};

cluster_matsubara_realspace_gf NormalStateClusterTransformer::transform_into_real_space(const cluster_matsubara_kspace_gf& G) const{
  cluster_matsubara_realspace_gf G_new(alps::gf::omega_r1_r2_sigma_gf(G.mesh1(), real_space_mesh(), real_space_mesh(), G.mesh3()));
  G_new.initialize();
  for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
    for (alps::gf::real_space_index i(0); i<G_new.mesh2().extent(); ++i) {
      for (alps::gf::real_space_index j(0); j<G_new.mesh3().extent(); ++j) {
        double phase = 0.;
        for (unsigned int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k(),d)*(cluster_coordinate(i(),d)-cluster_coordinate(j(),d));
        Backward::add(G_new, G, i(), j(), k(), cos(phase)/double(G.mesh2().extent()));
      }
    }
  }
  return G_new;
}
cluster_itime_realspace_gf NormalStateClusterTransformer::transform_into_real_space(const cluster_itime_kspace_gf& G) const{
  cluster_itime_realspace_gf G_new(alps::gf::itime_r1_r2_sigma_gf(G.mesh1(), real_space_mesh(), real_space_mesh(), G.mesh3()));
  G_new.initialize();
  for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
    for (alps::gf::real_space_index i(0); i<G_new.mesh2().extent(); ++i) {
      for (alps::gf::real_space_index j(0); j<G_new.mesh3().extent(); ++j) {
        double phase = 0.;
        for (int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k(),d)*(cluster_coordinate(i(),d)-cluster_coordinate(j(),d));
        Backward::add(G_new, G, i(), j(), k(), cos(phase)/double(G.mesh2().extent()));
      }
    }
  }
  return G_new;
}
cluster_matsubara_kspace_gf NormalStateClusterTransformer::transform_into_k_space(const cluster_matsubara_realspace_gf& G) const{
  cluster_matsubara_kspace_gf G_new(alps::gf::omega_k_sigma_gf(G.mesh1(), momentum_mesh(),G.mesh4()));
  G_new.initialize();
  for (alps::gf::momentum_index k(0); k<G_new.mesh2().extent(); ++k) {
    for (alps::gf::real_space_index i(0); i<G.mesh2().extent(); ++i) {
      for (alps::gf::real_space_index j(0); j<G.mesh3().extent(); ++j) {
        double phase = 0.;
        for (int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k(),d)*(cluster_coordinate(i(),d)-cluster_coordinate(j(),d));
        Forward::add(G_new, G, i(), j(), k(),  cos(phase)/double(G.mesh2().extent()));
      }
    }
  }
  return G_new;
}
cluster_itime_kspace_gf NormalStateClusterTransformer::transform_into_k_space(const cluster_itime_realspace_gf& G) const{
  cluster_itime_kspace_gf G_new(alps::gf::itime_k_sigma_gf(G.mesh1(), momentum_mesh(), G.mesh4()));
  G_new.initialize();
  for (alps::gf::momentum_index k(0); k<G_new.mesh2().extent(); ++k) {
    for (alps::gf::real_space_index i(0); i<G.mesh2().extent(); ++i) {
      for (alps::gf::real_space_index j(0); j<G.mesh3().extent(); ++j) {
        double phase = 0.;
        for (int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k(),d)*(cluster_coordinate(i(),d)-cluster_coordinate(j(),d));
        Forward::add(G_new, G, i(), j(), k(), cos(phase)/double(G.mesh2().extent()));
      }
    }
  }
  return G_new;
}

///Use Dyson to get the self-energy
void NormalStateClusterTransformer::dyson_get_self_energy(cluster_matsubara_kspace_gf &sigma, const cluster_matsubara_kspace_gf &G, const cluster_matsubara_kspace_gf &G0) const{
  for (alps::gf::matsubara_index w(0); w<G.mesh1().extent(); ++w) {
    for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
      for (alps::gf::index f(0); f<G.mesh3().extent(); ++f) {
        sigma(w,k,f)=1./G0(w,k,f)-1./G(w,k,f);
      }
    }
  }
}
///Use Dyson to get the G0
void NormalStateClusterTransformer::dyson_get_green0(cluster_matsubara_kspace_gf &G0, const cluster_matsubara_kspace_gf &sigma, const cluster_matsubara_kspace_gf &G) const{
  for (alps::gf::matsubara_index w(0); w<G.mesh1().extent(); ++w) {
    for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
      for (alps::gf::index f(0); f<G.mesh3().extent(); ++f) {
        G0(w,k,f)=1./(1./G(w,k,f)+sigma(w,k,f));
      }
    }
  }
}
///Use Dyson to get the G0
void NormalStateClusterTransformer::matsubara_k_sum(cluster_matsubara_kspace_gf &G, const cluster_matsubara_kspace_gf &sigma) const{
  for (alps::gf::matsubara_index w(0); w<G.mesh1().extent(); ++w) {
    for (alps::gf::momentum_index kc(0); kc<G.mesh2().extent(); ++kc) {
      for (alps::gf::index f(0); f<G.mesh3().extent(); ++f) {
        G(w,  kc, f) = 0.;
        std::complex<double> iomega_n(0, G.mesh1().points()[w()]);
        for (int kl=0; kl<n_lattice_momenta(); ++kl) {
          std::complex<double> zeta = iomega_n + mu_ + (f()%2 ? -h_ : h_);
          G(w, kc,f) += weight_kl_[kl]/(zeta-epsilon_kckl_[kc()][kl][f()]-sigma(w, kc,  f));
        }
        G(w, kc, f) /= n_cluster_cells();
      }
    }
  }
}

void NormalStateClusterTransformer::compute_integration_weight_and_dispersion(){
  std::vector<double> k(dimension());
  for (int kc=0; kc<n_cluster_sites(); ++kc) {
    for (int kl=0; kl<n_lattice_momenta(); ++kl) {
      for (unsigned int d=0; d<dimension(); ++d)
        k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
      for (unsigned int flavor=0; flavor<n_flavors_; ++flavor){
        epsilon_kckl_[kc][kl][flavor]=use_imbalance_?epsilon_flavor(k,flavor):use_phase_?epsilon_phi(k,phi_):epsilon(k);
      }
    }
  }
  for (int kl=0; kl<n_lattice_momenta(); ++kl) {
    weight_kl_[kl]= lattice_momentum(kl,dimension());
  }
}
void NormalStateClusterTransformer::compute_dispersion_weights(){
  epsav_             .resize(n_cluster_sites());
  epssqav_           .resize(n_cluster_sites());
  dispersion_moments_.resize(2);
  for (int kc=0; kc<n_cluster_sites(); ++kc) {
    for (int kl=0; kl<n_lattice_momenta(); ++kl) {
      std::vector<double> k(dimension());
      for (unsigned int d=0; d<dimension(); ++d)
        k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
      double weight = lattice_momentum(kl, dimension());
      double E = epsilon(k);
      epsav_  [kc]+= weight*E;
      epssqav_[kc]+= weight*E*E;
    }
    epsav_  [kc] /= n_cluster_cells();
    epssqav_[kc] /= n_cluster_cells();
    dispersion_moments_[0] += epsav_[kc];
    dispersion_moments_[1] += epssqav_[kc];
  }
  dispersion_moments_[0] /= n_cluster_sites();
  dispersion_moments_[1] /= n_cluster_sites();
}

double NormalStateClusterTransformer::epsilon(const std::vector<double> &k) const
{
  double eps = -half_filling_;
  for (unsigned int i=0; i<dimension(); ++i)
    eps -= 2*t_*cos(k[i]);
  if (t2_!=0.) {
    switch(dimension()) {
    case 1 :
      eps -= 2*t2_*cos(k[0]*2.);
      break;
    case 2 :
      eps -= 4*t2_*cos(k[0])*cos(k[1]);
      break;
    case 3 :
      eps -= 4*t2_*cos(k[0])*cos(k[1]);
      eps -= 4*t2_*cos(k[1])*cos(k[2]);
      eps -= 4*t2_*cos(k[0])*cos(k[2]);
      break;
    }
  }
  return eps;
}
void NormalStateClusterTransformer::setup_phase_vector(const alps::params &p){
  phi_.resize(dimension());
  for(int i=0;i<dimension();++i){
    std::stringstream phase_name; phase_name<<"dca.PHASE_"<<i;
    phi_[i]=p[phase_name.str()];
    phi_[i] = phi_[i]/periodicity(i);
  }
}

//dispersion, evaluated with flavor dependent t_
double NormalStateClusterTransformer::epsilon_flavor(const std::vector<double> &k, unsigned int flavor) const
{
  double eps = -half_filling_;

  if (flavor == 0){
    for (unsigned int i=0; i<dimension(); ++i)
      eps -= 2*t_*cos(k[i]);
  }
  if (flavor == 1){
    for (unsigned int i=0; i<dimension(); ++i)
      eps -= 2*t_*t_flavor_*cos(k[i]);
  }
  return eps;
}
double NormalStateClusterTransformer::epsilon_phi(const std::vector<double> &k, const std::vector<double> &phi) const{
  using namespace alps::numeric;
  return epsilon(k+phi);
}

