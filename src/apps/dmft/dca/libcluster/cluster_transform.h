/*****************************************************************************
 *
 * Copyright (C) 2007      by Sebastian Fuchs <fuchs@comp-phys.org>
 * Copyright (C) 2005-2026 by Emanuel Gull
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#pragma once

#include <Eigen/Dense>
#include <fstream>
#include "alps/gf/fourier.hpp"
#include "cluster_graph.h"
#include "cluster_helper.h"

#include "alps/gf/tail.hpp"
typedef alps::gf::omega_k_sigma_gf_with_tail cluster_matsubara_kspace_gf;
typedef alps::gf::itime_k_sigma_gf_with_tail cluster_itime_kspace_gf;
typedef alps::gf::omega_r1_r2_sigma_gf_with_tail cluster_matsubara_realspace_gf;
typedef alps::gf::itime_r1_r2_sigma_gf_with_tail cluster_itime_realspace_gf;

typedef alps::gf::omega_k_sigma1_sigma2_gf_with_tail cluster_nambu_matsubara_kspace_gf;
typedef alps::gf::itime_k_sigma1_sigma2_gf_with_tail cluster_nambu_itime_kspace_gf;
typedef alps::gf::omega_r1_r2_sigma1_sigma2_gf_with_tail cluster_nambu_matsubara_realspace_gf;
typedef alps::gf::itime_r1_r2_sigma1_sigma2_gf_with_tail cluster_nambu_itime_realspace_gf;

typedef alps::gf::omega_k_sigma1_sigma2_gf_with_tail double_cell_matsubara_kspace_gf;
typedef alps::gf::itime_k_sigma1_sigma2_gf_with_tail double_cell_itime_kspace_gf;


void define_dca_parameters(alps::params &p);

/// A table of values at cluster momentum kc plus lattice momentum kl.
class kckl_table{
public:
  kckl_table(int n_kc, int n_kl): n_kl_(n_kl), values_(std::size_t(n_kc)*n_kl, 0.){}
  double &operator()(int kc, int kl){ return values_[std::size_t(kc)*n_kl_+kl]; }
  double operator()(int kc, int kl) const{ return values_[std::size_t(kc)*n_kl_+kl]; }
private:
  int n_kl_;
  std::vector<double> values_;
};

//The Cluster transformer presents a framework for operator(), which does the Hilbert transform..
class ClusterTransformer : public ClusterHelper
{
public:
  ///Constructor: a cluster transformer has no internal components, everything is delegated to the ClusterHelper
  ClusterTransformer(const alps::params& params);
  ///virtual destructor
  virtual ~ClusterTransformer() {}

  ///mesh points for the k-space mesh of a DCA cluster transformer
  alps::gf::momentum_index_mesh momentum_mesh() const;
  ///mesh points for the real-space mesh of a DCA cluster transformer
  alps::gf::real_space_index_mesh real_space_mesh() const;


  ///dispersion relation for t and t' in 1, 2, and 3 d.
  double epsilon(const std::vector<double> &k) const;

  ///reading input file
  template<typename GF> bool read_initial_G0(GF &G0_omega) const;

protected:
  std::vector<double> epsav() const { return epsav_; };
  std::vector<double> epssqav() const { return epssqav_; };
  double dispersion_moments(const unsigned int i) const { return dispersion_moments_[i]; }
  ///get epsilon(k) and the simpson integration weights on the BZ
  void compute_integration_weight_and_dispersion();

  ///compute the center of gravity of the bands etc (used for high frequency analysis)
  void compute_dispersion_weights();

  ///the dispersion at cluster momentum kc plus lattice momentum kl; kl<0 means at the cluster momentum itself
  double dispersion(int kc, int kl) const;
  /// The moments m1=<H> and m2=<H^2> of the bare Hamiltonian H(kl) over the
  /// coarse-graining patch of cluster momentum kc, with the weights of the
  /// k-sum. The bare G0 = (<(iw - H - Sigma)^-1>^-1 + Sigma)^-1 then has the tail
  /// 1/iw + m1/(iw)^2 + m2/(iw)^3 for any Sigma. For the isolated cluster,
  /// H is taken at the cluster momentum only (hamiltonian(-1)).
  template<class M, class F> void bare_moments(int kc, bool isolated_cluster, F hamiltonian, M &m1, M &m2) const{
    if (isolated_cluster) {
      m1 = hamiltonian(-1);
      m2 = m1*m1;
      return;
    }
    m1.setZero(); m2.setZero();
    for (int kl=0; kl<n_lattice_momenta(); ++kl) {
      const double w = weight_kl_[kl]/n_cluster_cells();
      const M H = hamiltonian(kl);
      m1 += w*H;
      m2 += w*H*H;
    }
  }


  /// The starting G0 of the loop, and its Fourier transform: read from
  /// dca.G0_OMEGA_INPUT if that exists, otherwise computed by free (the
  /// noninteracting lattice) or isolated (the isolated cluster), as
  /// dca.INITIAL_G0 says, and then written to dca.G0_OMEGA_INPUT.
  template<typename GFW, typename GFT, typename FREE, typename ISOLATED>
  void initial_G0(GFW &G0_omega, GFT &G0_tau, FREE free, ISOLATED isolated) const{
    if(read_initial_G0(G0_omega)) return;
    if(initial_g0_=="FREE")
      free(G0_omega);
    else if (initial_g0_=="ISOLATED")
      isolated(G0_omega);
    else
      throw std::invalid_argument("choose FREE or ISOLATED as initial G0 parameters or provide an input hdf5 file");
    {
      alps::hdf5::archive ar(g0_input_path_, "w");
      G0_omega.save(ar, "/G0");
      std::ofstream G0_omega_input_file("G0_omega_input.dat"); G0_omega_input_file<<G0_omega;
    }
    alps::gf::fourier_frequency_to_time(G0_omega, G0_tau);
    std::ofstream G0_tau_input_file("G0_tau_input.dat"); G0_tau_input_file<<G0_tau;
  }
  /// Green's functions on the frequency and imaginary-time meshes of the loop,
  /// with the given index meshes after the momentum mesh.
  template<typename GFW, typename... INDEX> GFW make_matsubara_gf(INDEX... index) const{
    return GFW(typename GFW::gf_type(alps::gf::matsubara_positive_mesh(beta_, nfreq_), momentum_mesh(), index...));
  }
  template<typename GFT, typename... INDEX> GFT make_itime_gf(int n_tau, INDEX... index) const{
    return GFT(typename GFT::gf_type(alps::gf::itime_mesh(beta_, n_tau), momentum_mesh(), index...));
  }


  ///chemical potential
  double mu_;
  double half_filling_;
  ///t is the basic energy unit of the hopping. t2 is the tprime.
  double t_, t2_;
  double beta_;
  ///interaction is needed for high frequency information in the Fourier transform; same for the V
  double U_;
  double V_;
  int n_flavors_;
  int nfreq_;

  ///dispersion, binned by cluster momentum and lattice momentum (belonging to the cluster)
  kckl_table epsilon_kckl_;
  std::vector<double> weight_kl_;

  std::vector<double> epsav_;
  std::vector<double> epssqav_;
  std::vector<double> dispersion_moments_;
  std::string g0_input_path_;
  std::string initial_g0_;

};
class NormalStateClusterTransformer: public ClusterTransformer{
public:
  typedef cluster_matsubara_kspace_gf matsubara_gf;
  typedef cluster_itime_kspace_gf itime_gf;

  ///Constructor: a cluster transformer has no internal components, everything is delegated to the ClusterHelper
  NormalStateClusterTransformer(const alps::params& params): ClusterTransformer(params){}

  matsubara_gf make_matsubara_gf() const { return ClusterTransformer::make_matsubara_gf<matsubara_gf>(alps::gf::index_mesh(n_flavors_)); }
  itime_gf make_itime_gf(int n_tau) const { return ClusterTransformer::make_itime_gf<itime_gf>(n_tau, alps::gf::index_mesh(n_flavors_)); }
  void initial_G0(matsubara_gf &G0_omega, itime_gf &G0_tau) const{
    ClusterTransformer::initial_G0(G0_omega, G0_tau,
                                   [this](matsubara_gf &G){ compute_noninteracting_G0(G); },
                                   [this](matsubara_gf &G){ compute_isolated_cluster_G0(G); });
  }
  ///virtual destructor
  virtual ~NormalStateClusterTransformer() {}
  
  ///Use Dyson to get the self-energy
  void dyson_get_self_energy(cluster_matsubara_kspace_gf &sigma, const cluster_matsubara_kspace_gf &G, const cluster_matsubara_kspace_gf &G0) const;
  ///Use Dyson to get the G0
  void dyson_get_green0(cluster_matsubara_kspace_gf &G0, const cluster_matsubara_kspace_gf &sigma, const cluster_matsubara_kspace_gf &G) const;
  ///perform the k-sum to get G
  void matsubara_k_sum(cluster_matsubara_kspace_gf &G, const cluster_matsubara_kspace_gf &sigma) const;
  ///perform the k-sum to get G, adjust the chemical potential internally
  double matsubara_k_sum_adjust_mu(cluster_matsubara_kspace_gf &G, const cluster_matsubara_kspace_gf &sigma, double target_density);

  ///perform spin symmetrization of imaginary time k-space cluster Green's function
  void spin_symmetrize(cluster_itime_kspace_gf& G) const;
  ///perform spin symmetrization of matsubara k-space Green's function
  void spin_symmetrize(cluster_matsubara_kspace_gf& G) const;
  ///perform spin symmetrization of imaginary time k-space cluster Green's function
  void spin_symmetrize(cluster_itime_realspace_gf& G) const;
  ///perform spin symmetrization of matsubara k-space Green's function
  void spin_symmetrize(cluster_matsubara_realspace_gf& G) const;

  ///compute G0 from an atomic solution
  void compute_isolated_cluster_G0(cluster_matsubara_kspace_gf &G0_omega) const;
  ///compute G0 from a noninteracting solution
  void compute_noninteracting_G0(cluster_matsubara_kspace_gf &G0_omega) const;

  ///attach the DCA cluster high frequency behavior for the non-interacting GF
  void attach_bare_tail(cluster_matsubara_kspace_gf &G0, bool isolated_cluster=false) const;
  ///attach the DCA cluster high frequency behavior for the interacting GF
  void attach_interacting_tail(cluster_matsubara_kspace_gf &G, double U, const std::vector<double> &densities) const;

  //Real space and k-space transformers for the Green's functions
  cluster_matsubara_realspace_gf transform_into_real_space(const cluster_matsubara_kspace_gf& G) const;
  cluster_matsubara_kspace_gf transform_into_k_space(const cluster_matsubara_realspace_gf& G) const;
  cluster_itime_realspace_gf transform_into_real_space(const cluster_itime_kspace_gf& G) const;
  cluster_itime_kspace_gf transform_into_k_space(const cluster_itime_realspace_gf& G) const;


private:
  //symmetrization in momentum space - for the normal state GF
  template<class GF,class TAIL> void symmetrize_helper_k_space(alps::gf::three_index_gf_with_tail<GF,TAIL> &G,
                                                        const std::pair<std::set<pair_class_type>::const_iterator,
                                                        std::set<pair_class_type>::const_iterator> pair_classes) const;
  template<class GF, class TAIL> void symmetrize_helper_real_space(alps::gf::four_index_gf_with_tail<GF,TAIL> &G,
                                                        const std::pair<std::set<pair_class_type>::const_iterator,
                                                        std::set<pair_class_type>::const_iterator> pair_classes) const;
  ///given a chemical potential, compute the density (assuming a tail for the target density)
  double density_for_mu(double mu, const cluster_matsubara_kspace_gf &sigma, double target_density_per_site);

  /**
    *  optimized Fourier transformation from reciprocal space to real space
    *  The return Green's function is defined based on the GFT's first mesh type and value type.
    *  As the first the phase factor will be computed. Then using loop order based on the index order 
    *  the Fourier transform will performed. 
    *  
    *  @tparam GFT - type of K-space Green's to be transformed
    *
    *  @return - proper real-space Green's function for the GFT Green's function
    */
  template<typename GFT>
   alps::gf::four_index_gf_with_tail<alps::gf::four_index_gf<typename GFT::gf_type::value_type,
                  typename GFT::gf_type::mesh1_type,
                  alps::gf::real_space_index_mesh,
                  alps::gf::real_space_index_mesh,
                  alps::gf::index_mesh>,
                  alps::gf::three_index_gf<typename GFT::gf_type::value_type,
                                alps::gf::real_space_index_mesh,
                                alps::gf::real_space_index_mesh,
                                alps::gf::index_mesh> >
          transform_into_real_space_(const GFT & G) const{
    // init output Green's function
    alps::gf::four_index_gf_with_tail<alps::gf::four_index_gf<typename GFT::gf_type::value_type, // Green's function return value
                  typename GFT::gf_type::mesh1_type, // mesh type (frequency or time)
                  alps::gf::real_space_index_mesh,   // first real space index 
                  alps::gf::real_space_index_mesh,   // second real space index
                  alps::gf::index_mesh>,             // spin index
                  alps::gf::three_index_gf<typename GFT::gf_type::value_type, // tail specialization
                                alps::gf::real_space_index_mesh,
                                alps::gf::real_space_index_mesh,
                                alps::gf::index_mesh> >
       G_new(alps::gf::four_index_gf<typename GFT::gf_type::value_type,
                  typename GFT::gf_type::mesh1_type,
                  alps::gf::real_space_index_mesh,
                  alps::gf::real_space_index_mesh,
                  alps::gf::index_mesh> (G.mesh1(), real_space_mesh(), real_space_mesh(), G.mesh3()));
    G_new.initialize();
    // define Fourier phase factor as alps Green's function
    alps::gf::three_index_gf<double, alps::gf::momentum_index_mesh, alps::gf::real_space_index_mesh, alps::gf::real_space_index_mesh>
      phase(G.mesh2(), G_new.mesh2(), G_new.mesh3());
    // precompute fourier phase factor
    for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
      for (alps::gf::real_space_index i(0); i < G_new.mesh2().extent(); ++i) {
        for (alps::gf::real_space_index j(0); j < G_new.mesh3().extent(); ++j) {
          // accumulate phase factor into temporary variable
          double ph = 0.;
          for (unsigned int d=0; d<dimension(); ++d)
            ph += cluster_momentum(k(),d)*(cluster_coordinate(i(),d)-cluster_coordinate(j(),d));
          // store phase factor for <k, i, j> indices
          phase(k, i, j) = cos(ph) / (double)(G.mesh2().extent());
        }
      }
    }
    // iterate over indices according to the order in the Green's functions
    for (typename GFT::gf_type::mesh1_type::index_type w(0); w<G.mesh1().extent(); ++w) {
      for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
        for (alps::gf::real_space_index i(0); i<G_new.mesh2().extent(); ++i) {
          for (alps::gf::real_space_index j(0); j<G_new.mesh3().extent(); ++j) {
            for (alps::gf::index flavor(0); flavor<G.mesh3().extent(); ++flavor) {
              G_new(w, i, j, flavor) += phase(k, i, j) * G ( w, k, flavor);
            }
          }
        }
      }
    }
    // return new Green's function
    return G_new;
  }

public:
  template<class GF> void symmetrize_real_space(GF& G) const{return symmetrize_helper_real_space(G,pair_classes_real_space());}
  template<class GF> void symmetrize_k_space(GF& G) const{return symmetrize_helper_k_space(G,pair_classes_k_space());}
};

class DoubleCellClusterTransformer: public ClusterTransformer{
public:
  typedef double_cell_matsubara_kspace_gf matsubara_gf;
  typedef double_cell_itime_kspace_gf itime_gf;

  ///Constructor: a cluster transformer has no internal components, everything is delegated to the ClusterHelper
  DoubleCellClusterTransformer(const alps::params& params):
    ClusterTransformer(params),h_init_(params["dca.H_INIT"]), h_(params["dca.H"]),
    mu_stag_init_(params["dca.MU_STAG_INIT"]), mu_stag_(params["dca.MU_STAG"]), afm_symmetry_(params["dca.AFM_SELFCONSISTENCY"]), charge_order_symmetry_(params["dca.CO_SELFCONSISTENCY"]){
    if(afm_symmetry_ && mu_stag_!=0) throw std::invalid_argument("do not use staggered mu in AFM Self-consistency");
    if(charge_order_symmetry_ && h_!=0) throw std::invalid_argument("do not use magnetic field in CO Self-consistency");
    check_momentum_pairs();
  }
  ///virtual destructor
  virtual ~DoubleCellClusterTransformer() {}

  matsubara_gf make_matsubara_gf() const { return ClusterTransformer::make_matsubara_gf<matsubara_gf>(alps::gf::index_mesh(n_flavors_), alps::gf::index_mesh(n_flavors_)); }
  itime_gf make_itime_gf(int n_tau) const { return ClusterTransformer::make_itime_gf<itime_gf>(n_tau, alps::gf::index_mesh(n_flavors_), alps::gf::index_mesh(n_flavors_)); }
  void initial_G0(matsubara_gf &G0_omega, itime_gf &G0_tau) const{
    ClusterTransformer::initial_G0(G0_omega, G0_tau,
                                   [this](matsubara_gf &G){ compute_noninteracting_G0(G); },
                                   [this](matsubara_gf &G){ compute_isolated_cluster_G0(G); });
  }


  ///compute G0 from an atomic solution
  void compute_isolated_cluster_G0(double_cell_matsubara_kspace_gf &G0_omega) const;
  ///compute G0 from a noninteracting solution
  void compute_noninteracting_G0(double_cell_matsubara_kspace_gf &G0_omega) const;

  ///Use Dyson to get the self-energy
  void dyson_get_self_energy(double_cell_matsubara_kspace_gf &sigma, const double_cell_matsubara_kspace_gf &G, const double_cell_matsubara_kspace_gf &G0) const;
  ///Use Dyson to get the G0
  void dyson_get_green0(double_cell_matsubara_kspace_gf &G0, const double_cell_matsubara_kspace_gf &sigma, const double_cell_matsubara_kspace_gf &G) const;
  ///the k-sum
  void matsubara_k_sum(double_cell_matsubara_kspace_gf &G, const double_cell_matsubara_kspace_gf &sigma, double h, double mu_stag) const;
  ///the k-sum: short version, reroutes via long version by using h_ and mu_stag_
  void matsubara_k_sum(double_cell_matsubara_kspace_gf &G, const double_cell_matsubara_kspace_gf &sigma) const{
    matsubara_k_sum(G, sigma, h_, mu_stag_);
  }
  //the k-sum that will find the chem pot
  double matsubara_k_sum_adjust_mu(double_cell_matsubara_kspace_gf &G, const double_cell_matsubara_kspace_gf &sigma, double h, double mu_stag, double target_density);
  ///perform the k-sum to get G, adjust the chemical potential internally
  double matsubara_k_sum_adjust_mu(double_cell_matsubara_kspace_gf &G, const double_cell_matsubara_kspace_gf &sigma, double target_density){
    return matsubara_k_sum_adjust_mu(G, sigma, h_, mu_stag_,target_density);
  }

  ///perform spin symmetrization of matsubara k-space Green's function
  void spin_symmetrize(double_cell_matsubara_kspace_gf& G) const;
  ///average G over the classes of (K, K+Q) block elements that the cluster's point group relates
  void symmetrize_k_space(double_cell_matsubara_kspace_gf& G) const;

  ///take a cluster Green's function and transform it from k-space into real space
  cluster_itime_realspace_gf transform_into_real_space(const double_cell_itime_kspace_gf& G) const;

private:
  ///the blocks of G are (K, K+Q): check that momentum 2m+1 is momentum 2m plus Q
  void check_momentum_pairs() const;
  ///attach a tail to the noninteracting GF
  void attach_bare_tail(double_cell_matsubara_kspace_gf &G0, double h, double mu_stag, bool isolated_cluster=false) const;
  ///attach a tail to the interacting GF
  void attach_interacting_tail(double_cell_matsubara_kspace_gf &G, double U, double V, const std::vector<double> &densities) const;
  ///given a chemical potential, compute the density (assuming a tail for the target density)
  double density_for_mu(double mu, const double_cell_matsubara_kspace_gf &sigma, double target_density_per_site);

  ///initial magnetic field value
  const double h_init_;
  ///magnetic field
  const double h_;
  ///initial staggered chem pot value
  const double mu_stag_init_;
  ///staggered chem pot value
  const double mu_stag_;
  ///are we doing charge order symmetry?
  const bool charge_order_symmetry_;
  ///are we doing AFM symmertry?
  const bool afm_symmetry_;

};

class NambuClusterTransformer : public ClusterTransformer{
  ///Constructor: a cluster transformer has no internal components, everything is delegated to the ClusterHelper
public:
  NambuClusterTransformer(const alps::params& params);
  ///virtual destructor
  virtual ~NambuClusterTransformer(){}

  typedef cluster_nambu_matsubara_kspace_gf matsubara_gf;
  typedef cluster_nambu_itime_kspace_gf itime_gf;
  matsubara_gf make_matsubara_gf() const { return ClusterTransformer::make_matsubara_gf<matsubara_gf>(alps::gf::index_mesh(n_flavors_), alps::gf::index_mesh(n_flavors_)); }
  itime_gf make_itime_gf(int n_tau) const { return ClusterTransformer::make_itime_gf<itime_gf>(n_tau, alps::gf::index_mesh(n_flavors_), alps::gf::index_mesh(n_flavors_)); }
  void initial_G0(matsubara_gf &G0_omega, itime_gf &G0_tau) const{
    ClusterTransformer::initial_G0(G0_omega, G0_tau,
                                   [this](matsubara_gf &G){ compute_noninteracting_G0(G); },
                                   [this](matsubara_gf &G){ compute_isolated_cluster_G0(G); });
  }

  ///perform spin symmetrization of matsubara k-space Green's function
  void spin_symmetrize(cluster_nambu_matsubara_kspace_gf& G) const;

  ///compute G0 from an atomic solution
  void compute_isolated_cluster_G0(cluster_nambu_matsubara_kspace_gf &G0_omega) const;
  ///compute G0 from a noninteracting solution
  void compute_noninteracting_G0(cluster_nambu_matsubara_kspace_gf &G0_omega) const;

  ///Use Dyson to get the self-energy
  void dyson_get_self_energy(cluster_nambu_matsubara_kspace_gf &sigma, const cluster_nambu_matsubara_kspace_gf &G, const cluster_nambu_matsubara_kspace_gf &G0) const;
  ///Use Dyson to get the G0
  void dyson_get_green0(cluster_nambu_matsubara_kspace_gf &G0, const cluster_nambu_matsubara_kspace_gf &sigma, const cluster_nambu_matsubara_kspace_gf &G) const;
  ///perform the k-sum to get G
  void matsubara_k_sum(double_cell_matsubara_kspace_gf &G, const double_cell_matsubara_kspace_gf &sigma, double eta1, double eta2, double h) const;
  ///the k-sum: short version, reroutes via long version by using h_ and mu_stag_
  void matsubara_k_sum(double_cell_matsubara_kspace_gf &G, const double_cell_matsubara_kspace_gf &sigma) const{
    matsubara_k_sum(G, sigma, eta1_, eta2_, h_);
  }
  ///perform the k-sum to get G
  double matsubara_k_sum_adjust_mu(double_cell_matsubara_kspace_gf &G, const double_cell_matsubara_kspace_gf &sigma, double eta1, double eta2, double target_density);
  ///the k-sum: short version, reroutes via long version by using h_ and mu_stag_
  double matsubara_k_sum_adjust_mu(double_cell_matsubara_kspace_gf &G, const double_cell_matsubara_kspace_gf &sigma, double target_density) {
    return matsubara_k_sum_adjust_mu(G, sigma, eta1_, eta2_, target_density);
  }
  //Real space and k-space transformers for the Nambu Green's function
  cluster_nambu_matsubara_realspace_gf nambu_transform_into_real_space(const cluster_nambu_matsubara_kspace_gf& G) const;
  cluster_nambu_itime_realspace_gf nambu_transform_into_real_space(const cluster_nambu_itime_kspace_gf& G) const;

  void attach_bare_tail(cluster_nambu_matsubara_kspace_gf &G0, double h, double eta1, double eta2, bool isolated_cluster=false) const;
  ///the d-wave form factor cos kx - cos ky at cluster momentum kc plus lattice momentum kl; kl<0 means at the cluster momentum itself
  double dwave_prefactor(int kc, int kl) const;
  // TODO: h should also enter interacting tail
  void attach_interacting_tail(cluster_nambu_matsubara_kspace_gf &G, double U, const std::vector<double> &densities) const;
  void symmetrize_dwave_kspace(cluster_nambu_matsubara_kspace_gf &G) const;
  void symmetrize_k_space(cluster_nambu_matsubara_kspace_gf &G) const{return symmetrize_dwave_kspace(G);}

private:
  ///initial magnetic field value, Ferromagnet only
  const double h_init_;
  ///magnetic field, Ferromagnet only
  const double h_;

  ///superconducting bias field eta1, initial value
  double eta1_init_;
  ///superconducting bias field eta2, initial value
  double eta2_init_;
  ///superconducting bias field eta1
  double eta1_;
  ///superconducting bias field eta2
  double eta2_;
  ///buffered d-wave prefactor stored at the lattice momenta
  kckl_table dwave_prefactor_kckl_;
};

  

template<class GF,class TAIL> void NormalStateClusterTransformer::symmetrize_helper_k_space(alps::gf::three_index_gf_with_tail<GF,TAIL> &G,
                                                    const std::pair<std::set<pair_class_type>::const_iterator,
                                                    std::set<pair_class_type>::const_iterator> pair_classes) const
{
  std::set<pair_class_type>::const_iterator itc, itc_end;
  for(std::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    for (typename GF::mesh1_type::index_type freq(0); freq<G.mesh1().extent(); ++freq) {
      for (typename GF::mesh3_type::index_type flavor(0); flavor<G.mesh3().extent(); ++flavor) {
        typename GF::value_type av=0;
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          if(it->first != it->second) continue;
          av += G(freq, typename GF::mesh2_type::index_type(it->first), flavor);
        }
        av /= itc->size();
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          if(it->first != it->second) continue;
          G(freq, typename GF::mesh2_type::index_type(it->first), flavor) = av;
        }
      }
    }
  }
}
template<class GF,class TAIL> void NormalStateClusterTransformer::symmetrize_helper_real_space(alps::gf::four_index_gf_with_tail<GF,TAIL> &G,
                                                    const std::pair<std::set<pair_class_type>::const_iterator,
                                                    std::set<pair_class_type>::const_iterator> pair_classes) const
{
  std::set<pair_class_type>::const_iterator itc, itc_end;
  for(std::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    for (typename GF::mesh1_type::index_type freq(0); freq<G.mesh1().extent(); ++freq) {
      for (typename GF::mesh4_type::index_type flavor(0); flavor<G.mesh4().extent(); ++flavor) {
        typename GF::value_type av=0;
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          av += G(freq, typename GF::mesh2_type::index_type(it->first),  typename GF::mesh3_type::index_type(it->second), flavor);
        }
        av /= itc->size();
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          G(freq, typename GF::mesh2_type::index_type(it->first), typename GF::mesh3_type::index_type(it->second), flavor) = av;
        }
      }
    }
  }
}

///lightweight class for operating 2x2 matrices
class two_by_two{
public:
  ///construct with four values
  two_by_two(std::complex<double> v11,std::complex<double> v12, std::complex<double> v21, std::complex<double> v22){ val[0]=v11; val[1]=v12; val[2]=v21; val[3]=v22;}
  ///access value
  const std::complex<double> &operator()(int p,int q) const{ return val[2*p+q];}
  std::complex<double> &operator()(int p,int q) { return val[2*p+q];}
  ///invert the entire matrix
  void invert(){
    std::complex<double> det_inv=1./(val[0]*val[3]-val[2]*val[1]);
    std::complex<double> tmp=val[0];
    val[0]=val[3]*det_inv;
    val[3]=tmp*det_inv;
    val[1]=-val[1]*det_inv;
    val[2]=-val[2]*det_inv;
  }
  two_by_two operator*(const two_by_two &rhs) const{
    return two_by_two(val[0]*rhs.val[0]+val[1]*rhs.val[2], val[0]*rhs.val[1]+val[1]*rhs.val[3], val[2]*rhs.val[0]+val[3]*rhs.val[2], val[2]*rhs.val[1]+val[3]*rhs.val[3]);
  }
  two_by_two operator-(const two_by_two &rhs) const{
    return two_by_two(val[0]-rhs.val[0],val[1]-rhs.val[1], val[2]-rhs.val[2],val[3]-rhs.val[3]);
  }
  two_by_two operator*(const std::complex<double> &v) const{
    return two_by_two(val[0]*v,val[1]*v, val[2]*v,val[3]*v);
  }

private:
  std::complex<double> val[4];
};
