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
  p.define<double>("t_flavor", 0, "Hopping for mass imbalance for second orbital");
  p.define<bool>("dca.IMBALANCE", false, "set to true for mass imbalance simulation");
  p.define<bool>("dca.PHASE", false, "set to true for simulation with phase");
  p.define<bool>("sc.NAMBU", false, "set to true to simulate in nambu space.");
  p.define<double>("dca.PHASE_0", 0, "Phase in direction 0 ");
  p.define<double>("dca.PHASE_1", 0, "Phase in direction 1 ");
  p.define<double>("dca.PHASE_2", 0, "Phase in direction 2 ");
  p.define<double>("dca.ETA1_INIT", 0, "d-wave initial bias field eta1");
  p.define<double>("dca.ETA2_INIT", 0, "d-wave initial bias field eta2");
  p.define<double>("dca.ETA1", 0, "d-wave bias field eta1");
  p.define<double>("dca.ETA2", 0, "d-wave bias field eta2");
  p.define<std::string>("dca.G0_OMEGA_INPUT", "Location of the input G0 file");
  p.define<std::string>("dca.G0_OMEGA_OUTPUT", "Location of the output G0 file for continuing iterations");
  p.define<std::string>("dca.INITIAL_G0", "FREE", "choose ISOLATED or FREE for initial G0s");
}
ClusterTransformer::ClusterTransformer(const alps::params& params) : ClusterHelper(params),
    mu_(params["MU"]),
    h_(0),
    half_filling_(0.),
    t_(params["t"]),
    t2_(params["tprime"]),
    t_flavor_(params["t_flavor"]),
    n_flavors_(params["FLAVORS"]),
    use_phase_(params["dca.PHASE"]),
    use_imbalance_(params["dca.IMBALANCE"]),
    epsilon_kckl_(boost::extents[n_cluster_sites()][n_lattice_momenta()][n_flavors_]),
    weight_kl_(n_lattice_momenta()),
    nfreq_(params["NMATSUBARA"]),
    beta_(params["BETA"]),
    g0_input_path_(params["dca.G0_OMEGA_INPUT"].as<std::string>()),
    initial_g0_(params["dca.INITIAL_G0"].as<std::string>())
{
  compute_integration_weight_and_dispersion();
  compute_dispersion_weights();
  if(use_phase_) setup_phase_vector(params);
}
NambuClusterTransformer::NambuClusterTransformer(const alps::params& params):ClusterTransformer(params),
    eta1_init_(params["dca.ETA1_INIT"]),
    eta2_init_(params["dca.ETA2_INIT"]),
    eta1_(params["dca.ETA1"]),
    eta2_(params["dca.ETA2"]),
    dwave_prefactor_kckl_(boost::extents[n_cluster_sites()][n_lattice_momenta()])
{
  if(dimension()!=2) throw std::logic_error("careful. self consistency has d=2 symmetry factors built in! do you know what you're doing?");

  for (int kc=0; kc<n_cluster_sites(); ++kc) {
    for (int kl=0; kl<n_lattice_momenta(); ++kl) {
      std::vector<double> k(dimension());
      for (unsigned int d=0; d<dimension(); ++d)
        k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
      dwave_prefactor_kckl_[kc][kl]= cos(k[0]) - cos(k[1]); //cos kx - cos ky
    }
  }
}

template<typename GF> bool ClusterTransformer::read_initial_G0(GF &G0_omega) const{
  std::cout<<"Reading initial G0_omega from: "<<g0_input_path_<<std::endl;
  {
    std::ifstream in_file(g0_input_path_.c_str());
    if(!in_file.good()){
      std::cout<<"WARNING: Input file could not be found, computing initial solution."<<std::endl;
      return false;
    }
  }

  alps::hdf5::archive ar(g0_input_path_);
  G0_omega.load(ar, "/G0");
  return true;
}
bool NormalStateClusterTransformer::read_initial_G0(cluster_matsubara_kspace_gf &G0_omega) const{
  return ClusterTransformer::read_initial_G0(G0_omega);
}
bool NambuClusterTransformer::read_initial_G0(cluster_nambu_matsubara_kspace_gf &G0_omega) const{
  return ClusterTransformer::read_initial_G0(G0_omega);
}

void NormalStateClusterTransformer::compute_isolated_cluster_G0(cluster_matsubara_kspace_gf &G0_omega) const{
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
  attach_bare_tail(G0_omega);
}
void NambuClusterTransformer::compute_isolated_cluster_G0(cluster_nambu_matsubara_kspace_gf &G0_omega) const{
  std::vector<double> k(dimension());

  for (alps::gf::momentum_index kc(0); kc<n_cluster_sites(); ++kc) {
    for (alps::gf::matsubara_index freq(0); freq<nfreq_; ++freq) {
      std::complex<double> iomega_n(0, G0_omega.mesh1().points()[freq()]);

      for (unsigned int d=0; d<dimension(); ++d)
        k[d] = cluster_momentum(kc(),d);
      double dwave_prefactor= cos(k[0]) - cos(k[1]); //cos kx - cos ky
      double prefactor=dwave_prefactor;

      //note the pairing field eta + i eta2 on the 0,1 component, eta - i eta2 on the 1,0 component.
      std::complex<double> G0_omega_latt_inv[4]=
      {iomega_n + mu_ -epsilon(k),
          - prefactor*std::complex<double>( eta1_init_, eta2_init_),
          - prefactor*std::complex<double>( eta1_init_,-eta2_init_),
          iomega_n - mu_ +epsilon(k)};
      std::complex<double> G0omega_latt_inv_det=G0_omega_latt_inv[0]*G0_omega_latt_inv[3]-G0_omega_latt_inv[1]*G0_omega_latt_inv[2];
      G0_omega(freq, kc, alps::gf::index(0),alps::gf::index(0)) =  G0_omega_latt_inv[3]/G0omega_latt_inv_det;
      G0_omega(freq, kc, alps::gf::index(0),alps::gf::index(1)) =- G0_omega_latt_inv[1]/G0omega_latt_inv_det;
      G0_omega(freq, kc, alps::gf::index(1),alps::gf::index(0)) =- G0_omega_latt_inv[2]/G0omega_latt_inv_det;
      G0_omega(freq, kc, alps::gf::index(1),alps::gf::index(1)) =  G0_omega_latt_inv[0]/G0omega_latt_inv_det;
    }
  }
  attach_bare_tail(G0_omega);
}
void NormalStateClusterTransformer::compute_noninteracting_G0(cluster_matsubara_kspace_gf &G0_omega) const{
  cluster_matsubara_kspace_gf zero_selfenergy(G0_omega);
  zero_selfenergy.initialize();
  matsubara_k_sum(G0_omega, zero_selfenergy);
  attach_bare_tail(G0_omega);
}
void NambuClusterTransformer::compute_noninteracting_G0(cluster_nambu_matsubara_kspace_gf &G0_omega) const{
  cluster_nambu_matsubara_kspace_gf zero_selfenergy(G0_omega);
  zero_selfenergy.initialize();
  matsubara_k_sum(G0_omega, zero_selfenergy);
  attach_bare_tail(G0_omega);
}

template<typename GFW, typename GFT> void ClusterTransformer::initial_G0(GFW &G0_omega, GFT &G0_tau)const{
  if(read_initial_G0(G0_omega)) return;

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
    alps::gf::fourier_frequency_to_time(G0_omega, G0_tau);
    std::ofstream G0_tau_input_file("G0_tau_input.dat"); G0_tau_input_file<<G0_tau;
  }
}

void NormalStateClusterTransformer::attach_bare_tail(cluster_matsubara_kspace_gf &G0) const{
  alps::gf::k_sigma_gf c1(G0.mesh2(), G0.mesh3());
  alps::gf::k_sigma_gf c2(G0.mesh2(), G0.mesh3());
  alps::gf::k_sigma_gf c3(G0.mesh2(), G0.mesh3());

  for(alps::gf::index f(0); f<n_flavors_; f++) {
    for(alps::gf::momentum_index k(0); k<n_cluster_sites(); ++k) {
      c1(k,f) = 1.;
      c2(k,f) = epsav_[k()] - mu_;
      c3(k,f) = epssqav_[k()] - 2*mu_*epsav_[k()] + mu_*mu_;
    }
  }
  G0.set_tail(1,c1);
  G0.set_tail(2,c2);
  G0.set_tail(3,c3);
}
void NormalStateClusterTransformer::attach_interacting_tail(cluster_matsubara_kspace_gf &G, double U, const std::vector<double> &densities) const{
  alps::gf::k_sigma_gf c1(G.mesh2(), G.mesh3());
  alps::gf::k_sigma_gf c2(G.mesh2(), G.mesh3());
  alps::gf::k_sigma_gf c3(G.mesh2(), G.mesh3());
  double density=(densities[0]+densities[1])/2.;

  for(alps::gf::index f(0); f<n_flavors_; f++) {
    for(alps::gf::momentum_index k(0); k<n_cluster_sites(); ++k) {
      c1(k,f) = 1.;
      c2(k,f) = epsav_[k()] - mu_+ U*density;
      c3(k,f) = epssqav_[k()] - 2*mu_*epsav_[k()] + mu_*mu_+ 2.*U*density*(epsav_[k()]-mu_) + U*U*density;
    }
  }
  G.set_tail(1,c1);
  G.set_tail(2,c2);
  G.set_tail(3,c3);
}

alps::gf::momentum_index_mesh ClusterTransformer::momentum_mesh() const{
  alps::gf::momentum_index_mesh::container_type points(boost::extents[n_cluster_sites()][dimension()]);
  for(int i=0;i<n_cluster_sites();++i){
    for(int j=0;j<dimension();++j){
      points[i][j]=cluster_momentum(i,j);
    }
  }
  return alps::gf::momentum_index_mesh(points);
}
alps::gf::real_space_index_mesh ClusterTransformer::real_space_mesh() const{
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
///perform spin symmetrization of matsubara k-space Green's function in Nambu space
void NambuClusterTransformer::spin_symmetrize(cluster_nambu_matsubara_kspace_gf& G) const{
  alps::gf::index up(0);
  alps::gf::index dn(1);
  if(G.mesh3().extent()!=2) throw std::runtime_error("mesh should only have two flavors.");
  for (alps::gf::matsubara_index w=alps::gf::matsubara_index(0); w<G.mesh1().extent(); ++w) {
    for (alps::gf::momentum_index k=alps::gf::momentum_index(0); k<G.mesh2().extent(); ++k) {
      //symmetrize the diagonal G
      G(w, k, up,up)=std::complex<double>(0.5*(G(w,k,up,up).real()-G(w,k,dn,dn).real()),0.5*(G(w,k,up,up).imag()+G(w,k,dn,dn).imag()));
      G(w, k, dn,dn)=std::complex<double>(-G(w,k,up,up).real(),G(w,k,dn,dn).imag());
      //symmetrize the  F.
      G(w, k, up, dn)= (G(w, k, up, dn)+std::conj(G(w, k, dn, up)))/2.;
      G(w, k, dn, up)= std::conj(G(w, k, up, dn));
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

///TODO: we could merge at least two of these, potentially four of these into one k-space FT
cluster_matsubara_realspace_gf NormalStateClusterTransformer::transform_into_real_space(const cluster_matsubara_kspace_gf& G) const{
  cluster_matsubara_realspace_gf G_new(alps::gf::omega_r1_r2_sigma_gf(G.mesh1(), real_space_mesh(), real_space_mesh(), G.mesh3()));
  G_new.initialize();
  for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
    for (alps::gf::real_space_index i(0); i<G_new.mesh2().extent(); ++i) {
      for (alps::gf::real_space_index j(0); j<G_new.mesh3().extent(); ++j) {
        double phase = 0.;
        for (unsigned int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k(),d)*(cluster_coordinate(i(),d)-cluster_coordinate(j(),d));
        for (alps::gf::matsubara_index w(0); w<G.mesh1().extent(); ++w) {
          for (alps::gf::index flavor(0); flavor<G.mesh3().extent(); ++flavor) {
            G_new(w, i, j, flavor) += (cos(phase)/(double)(G.mesh2().extent()))* G(w, k, flavor);
          }
        }
      }
    }
  }
  return G_new;
}
cluster_nambu_matsubara_realspace_gf NambuClusterTransformer::transform_into_real_space(const cluster_nambu_matsubara_kspace_gf& G) const{
  cluster_nambu_matsubara_realspace_gf G_new(alps::gf::omega_r1_r2_sigma1_sigma2_gf(G.mesh1(), real_space_mesh(), real_space_mesh(), G.mesh3(), G.mesh4()));
  G_new.initialize();
  for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
    for (alps::gf::real_space_index i(0); i<G_new.mesh2().extent(); ++i) {
      for (alps::gf::real_space_index j(0); j<G_new.mesh3().extent(); ++j) {
        double phase = 0.;
        for (unsigned int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k(),d)*(cluster_coordinate(i(),d)-cluster_coordinate(j(),d));
        for (alps::gf::matsubara_index w(0); w<G.mesh1().extent(); ++w) {
        for (alps::gf::index flavor1(0); flavor1<G.mesh3().extent(); ++flavor1) {
          for (alps::gf::index flavor2(0); flavor2<G.mesh3().extent(); ++flavor2) {
            G_new(w, i, j, flavor1, flavor2) += (cos(phase)/(double)(G.mesh2().extent()))* G(w, k, flavor1, flavor2);
          }
        }
        }
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
        for (alps::gf::itime_index t(0); t<G.mesh1().extent(); ++t) {
        for (alps::gf::index flavor(0); flavor<G.mesh3().extent(); ++flavor) {
          G_new(t, i, j, flavor) += (cos(phase)/(double)(G.mesh2().extent()))* G(t, k, flavor);
        }
        }
      }
    }
  }
  return G_new;
}
cluster_nambu_itime_realspace_gf NambuClusterTransformer::transform_into_real_space(const cluster_nambu_itime_kspace_gf& G) const{
  cluster_nambu_itime_realspace_gf G_new(alps::gf::itime_r1_r2_sigma1_sigma2_gf(G.mesh1(), real_space_mesh(), real_space_mesh(), G.mesh3(), G.mesh4()));
  G_new.initialize();
  for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
    for (alps::gf::real_space_index i(0); i<G_new.mesh2().extent(); ++i) {
      for (alps::gf::real_space_index j(0); j<G_new.mesh3().extent(); ++j) {
        double phase = 0.;
        for (int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k(),d)*(cluster_coordinate(i(),d)-cluster_coordinate(j(),d));
        for (alps::gf::itime_index t(0); t<G.mesh1().extent(); ++t) {
        for (alps::gf::index flavor1(0); flavor1<G.mesh3().extent(); ++flavor1) {
          for (alps::gf::index flavor2(0); flavor2<G.mesh3().extent(); ++flavor2) {
            G_new(t, i, j, flavor1, flavor2) += (cos(phase)/(double)(G.mesh2().extent()))* G(t, k, flavor1, flavor2);
          }
        }
        }
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

  //Set high frequency tails for self energy
  alps::gf::k_sigma_gf s0(G.mesh2(), G.mesh3());
  alps::gf::k_sigma_gf s1(G.mesh2(), G.mesh3());

  //Get the G tail coeffs
  alps::gf::k_sigma_gf c2(G.mesh2(), G.mesh3());
  alps::gf::k_sigma_gf c3(G.mesh2(), G.mesh3());
  c2 = G.tail(2);
  c3 = G.tail(3);

  //Get the G0 tail coeffs
  alps::gf::k_sigma_gf c02(G.mesh2(), G.mesh3());
  alps::gf::k_sigma_gf c03(G.mesh2(), G.mesh3());
  c02 = G0.tail(2);
  c03 = G0.tail(3);

  for(alps::gf::index f(0); f<n_flavors_; f++) {
    for(alps::gf::momentum_index k(0); k<n_cluster_sites(); ++k) {
      s0(k,f) = c2(k,f) - c02(k,f);
      s1(k,f) = c3(k,f) - c03(k,f) + c02(k,f)*c02(k,f) - c2(k,f)*c2(k,f);
    }
  }
  sigma.set_tail(0,s0);
  sigma.set_tail(1,s1);
}
///Use Dyson to get the self-energy in Nambu space
void NambuClusterTransformer::dyson_get_self_energy(cluster_nambu_matsubara_kspace_gf &sigma, const cluster_nambu_matsubara_kspace_gf &G, const cluster_nambu_matsubara_kspace_gf &G0) const{
  const alps::gf::index up(0);
  const alps::gf::index dn(1);
  for (alps::gf::matsubara_index w(0); w<G.mesh1().extent(); ++w) {
    for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
      std::complex<double> g0_block[4]={G0(w, k, up,up),G0(w, k, up,dn),G0(w, k, dn,up),G0(w, k, dn,dn)};
      std::complex<double>  g_block[4]={ G(w, k, up,up), G(w, k, up,dn), G(w, k, dn,up), G(w, k, dn,dn)};
      std::complex<double>      det_g0=g0_block[0]*g0_block[3]-g0_block[1]*g0_block[2];
      std::complex<double>      det_g = g_block[0]* g_block[3]- g_block[1]* g_block[2];

      sigma(w,k,up,up)= g0_block[3]/det_g0 - g_block[3]/det_g;
      sigma(w,k,up,dn)=-g0_block[1]/det_g0 + g_block[1]/det_g;
      sigma(w,k,dn,up)=-g0_block[2]/det_g0 + g_block[2]/det_g;
      sigma(w,k,dn,dn)= g0_block[0]/det_g0 - g_block[0]/det_g;
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
  alps::gf::k_sigma_gf tail_1(G0.mesh2(), G0.mesh3());
  alps::gf::k_sigma_gf tail_2(G0.mesh2(), G0.mesh3());
  for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
    for (alps::gf::index f(0); f<G.mesh3().extent(); ++f) {
      tail_1(k,f)=1.;
      tail_2(k,f)=dispersion_moments(k())-mu_;
    }
  }
  G0.set_tail(1, tail_1);
  G0.set_tail(2, tail_2);
}
///Use Dyson to get the G0 in Nambu space
void NambuClusterTransformer::dyson_get_green0(cluster_nambu_matsubara_kspace_gf &G0, const cluster_nambu_matsubara_kspace_gf &sigma, const cluster_nambu_matsubara_kspace_gf &G) const{
  const alps::gf::index up(0);
  const alps::gf::index dn(1);
  for (alps::gf::matsubara_index w(0); w<G.mesh1().extent(); ++w) {
    for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
      std::complex<double> glatt_block[4]={ G(w, k, up,up), G(w, k, up, dn), G(w, k, dn, up), G(w, k, dn, dn)};
      std::complex<double> det_glatt     = glatt_block[0]* glatt_block[3]- glatt_block[1]* glatt_block[2];
      std::complex<double> glatt_inv[4]  ={glatt_block[3]/det_glatt,-glatt_block[1]/det_glatt,-glatt_block[2]/det_glatt,glatt_block[0]/det_glatt};
      std::complex<double> g0_inv[4]     ={glatt_inv[0]+sigma(w,k,up,up), glatt_inv[1]+sigma(w,k,up,dn),glatt_inv[2]+sigma(w,k,dn,up),glatt_inv[3]+sigma(w,k,dn,dn)};
      std::complex<double> det_g0inv     =g0_inv[0]*g0_inv[3]-g0_inv[1]*g0_inv[2];

      G0(w,k,up,up) = g0_inv[3]/det_g0inv;
      G0(w,k,up,dn) =-g0_inv[1]/det_g0inv;
      G0(w,k,dn,up) =-g0_inv[2]/det_g0inv;
      G0(w,k,dn,dn) = g0_inv[0]/det_g0inv;
    }
  }
}

///Compute the DMFT/DCA k-summation
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
///Compute the DMFT/DCA k-summation
void NambuClusterTransformer::matsubara_k_sum(cluster_nambu_matsubara_kspace_gf &G, const cluster_nambu_matsubara_kspace_gf &sigma) const{
  alps::gf::index up(0);
  alps::gf::index dn(1);
  for (alps::gf::momentum_index kc(0); kc<G.mesh2().extent(); ++kc) {
    for (alps::gf::matsubara_index w(0); w<G.mesh1().extent(); ++w) {
      G(w, kc, up,up) = 0.;
      G(w, kc, up,dn) = 0.;
      G(w, kc, dn,up) = 0.;
      G(w, kc, dn,dn) = 0.;
      std::complex<double> iomega_n(0, (2.*w()+1.)*M_PI/beta_);
      for (int kl=0; kl<n_lattice_momenta(); ++kl) {
        double prefactor=dwave_prefactor_kckl_[kc()][kl];
        //this is the self consistency formula!!!
        //note the pairing field eta + i eta2 on the 0,1 component, eta - i eta2 on the 1,0 component.
        std::complex<double> G_omega_latt_inv[4]=
        {iomega_n + mu_ -epsilon_kckl_[kc()][kl][0]-sigma(w, kc, up,up),
            -sigma(w, kc, up,dn) - prefactor*std::complex<double>( eta1_, eta2_),
            -sigma(w, kc, dn,up) - prefactor*std::complex<double>( eta1_,-eta2_),
            iomega_n - mu_ +epsilon_kckl_[kc()][kl][0]-sigma(w, kc, dn,dn)};

        std::complex<double> G0omega_latt_inv_det=G_omega_latt_inv[0]*G_omega_latt_inv[3]-G_omega_latt_inv[1]*G_omega_latt_inv[2];
        G(w, kc, up,up) += weight_kl_[kl]*G_omega_latt_inv[3]/G0omega_latt_inv_det;
        G(w, kc, up,dn) -= weight_kl_[kl]*G_omega_latt_inv[1]/G0omega_latt_inv_det;
        G(w, kc, dn,up) -= weight_kl_[kl]*G_omega_latt_inv[2]/G0omega_latt_inv_det;
        G(w, kc, dn,dn) += weight_kl_[kl]*G_omega_latt_inv[0]/G0omega_latt_inv_det;
      }
      G(w, kc, up,up) /= n_cluster_cells();
      G(w, kc, up,dn) /= n_cluster_cells();
      G(w, kc, dn,up) /= n_cluster_cells();
      G(w, kc, dn,dn) /= n_cluster_cells();
    }
  }
}

void ClusterTransformer::compute_integration_weight_and_dispersion(){
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
void ClusterTransformer::compute_dispersion_weights(){
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

double ClusterTransformer::epsilon(const std::vector<double> &k) const
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
void ClusterTransformer::setup_phase_vector(const alps::params &p){
  phi_.resize(dimension());
  for(int i=0;i<dimension();++i){
    std::stringstream phase_name; phase_name<<"dca.PHASE_"<<i;
    phi_[i]=p[phase_name.str()];
    phi_[i] = phi_[i]/periodicity(i);
  }
}

//dispersion, evaluated with flavor dependent t_
double ClusterTransformer::epsilon_flavor(const std::vector<double> &k, unsigned int flavor) const
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
double ClusterTransformer::epsilon_phi(const std::vector<double> &k, const std::vector<double> &phi) const{
  using namespace alps::numeric;
  return epsilon(k+phi);
}
///merge this with normal state one?
void NambuClusterTransformer::attach_bare_tail(cluster_nambu_matsubara_kspace_gf &G0) const{
  alps::gf::k_sigma1_sigma2_gf c1(G0.mesh2(), G0.mesh3(), G0.mesh4());
  alps::gf::k_sigma1_sigma2_gf c2(G0.mesh2(), G0.mesh3(), G0.mesh4());
  alps::gf::k_sigma1_sigma2_gf c3(G0.mesh2(), G0.mesh3(), G0.mesh4());

  for(alps::gf::index f1(0); f1<n_flavors_; f1++) {
    for(alps::gf::index f2(0); f2<n_flavors_; f2++) {
      for(alps::gf::momentum_index k(0); k<n_cluster_sites(); ++k) {
        c1(k,f1,f2) = f1()==f2()?1.:0.;
        c2(k,f1,f2) = (f1()==0 && f2()==0)? epsav_[k()] - mu_: (f1()==1 && f2()==1)? -(epsav_[k()] - mu_):0.;
        c3(k,f1,f2) = 0.;
      }
    }
  }
  G0.set_tail(1,c1);
  G0.set_tail(2,c2);
  G0.set_tail(3,c3);
}
///merge this with normal state one?
void NambuClusterTransformer::attach_interacting_tail(cluster_nambu_matsubara_kspace_gf &G, double U, const std::vector<double> &densities) const{
  alps::gf::k_sigma1_sigma2_gf c1(G.mesh2(), G.mesh3(), G.mesh4());
  alps::gf::k_sigma1_sigma2_gf c2(G.mesh2(), G.mesh3(), G.mesh4());
  alps::gf::k_sigma1_sigma2_gf c3(G.mesh2(), G.mesh3(), G.mesh4());
  double density=(densities[0]+densities[1])/2.;

  for(alps::gf::index f1(0); f1<n_flavors_; f1++) {
    for(alps::gf::index f2(0); f2<n_flavors_; f2++) {
      for(alps::gf::momentum_index k(0); k<n_cluster_sites(); ++k) {
        c1(k,f1,f2) = f1()==f2()?1.:0.;
        c2(k,f1,f2) = (f1()==0 && f2()==0)? epsav_[k()] - mu_+ U*density: (f1()==1 && f2()==1)? -(epsav_[k()] - mu_ + U*density):0.;
        c3(k,f1,f2) = 0.;
      }
    }
  }
  G.set_tail(1,c1);
  G.set_tail(2,c2);
  G.set_tail(3,c3);
}

void NambuClusterTransformer::symmetrize_dwave_kspace(cluster_nambu_matsubara_kspace_gf &G) const
{
  //symmetrize the diagonal elements of G
  std::set<pair_class_type>::const_iterator itc, itc_end;
  for(boost::tie(itc, itc_end)=pair_classes_k_space(); itc!=itc_end; itc++) {
    for (alps::gf::index f(0); f<2; ++f) {
      for (alps::gf::matsubara_index w(0); w<nfreq_; ++w) {
        std::complex<double> av=0;
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          if(it->first != it-> second) continue;
          av += G(w, alps::gf::momentum_index(it->first), f,f);
        }
        av /= itc->size();
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          if(it->first != it-> second) continue;
          G(w, alps::gf::momentum_index(it->first), f, f) = av;
        }
      }
    }
  }

  //symmetrize the off-diagonal elements of G. Note that we assume that a the lattice symmetry can further be broken by a d-wave term
  std::vector<double> dwave_prefactor(G.mesh2().extent());
  for(unsigned int i=0;i<G.mesh2().extent();++i){
    double dwave=cos(cluster_momentum(i,0))-cos(cluster_momentum(i,1)); // <- this is the d-wave term
    if(dwave < -1.e-5) dwave_prefactor[i]=-1;
    else if(dwave >  1.e-5) dwave_prefactor[i]=1;
    else dwave_prefactor[i]=0.;
  }
  for(boost::tie(itc, itc_end)=pair_classes_k_space(); itc!=itc_end; itc++) {
    //std::cout<<"new symmetry: "<<std::endl;
    for (alps::gf::matsubara_index w(0); w<nfreq_; ++w) {
      std::complex<double> av01=0;
      std::complex<double> av10=0;
      for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
        if(it->first != it-> second) continue;
        av01+= dwave_prefactor[it->first]*G(w, alps::gf::momentum_index(it->first), alps::gf::index(0),alps::gf::index(1));
        av10+= dwave_prefactor[it->first]*G(w, alps::gf::momentum_index(it->first), alps::gf::index(1),alps::gf::index(0));
      }
      av01/= itc->size();
      av10/= itc->size();
      //if(freq==0){ std::cout<<"setting average to: "<<av01<<" "<<av10<<std::endl; }
      for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
        if(it->first != it-> second) continue;
        G(w, alps::gf::momentum_index(it->first), alps::gf::index(0),alps::gf::index(1)) = dwave_prefactor[it->first]*av01;
        G(w, alps::gf::momentum_index(it->first), alps::gf::index(1),alps::gf::index(0)) = dwave_prefactor[it->first]*av10;
      }
    }
  }
}

//create some instances for libcluster:
template void ClusterTransformer::initial_G0(cluster_nambu_matsubara_kspace_gf &G0_omega, cluster_nambu_itime_kspace_gf &G0_tau)const;
template void ClusterTransformer::initial_G0(cluster_matsubara_kspace_gf &G0_omega, cluster_itime_kspace_gf &G0_tau)const;
