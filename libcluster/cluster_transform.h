/*****************************************************************************
 *
 * Copyright (C) 2007      by Sebastian Fuchs <fuchs@comp-phys.org>
 * Copyright (C) 2007-2016 by Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/

#pragma once

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

  template<typename GFW, typename GFT> void initial_G0(GFW &G0_omega, GFT &G0_tau) const;

  ///dispersion relation for t and t' in 1, 2, and 3 d.
  double epsilon(const std::vector<double> &k) const;
  ///mass imbalance spin-dependent dispersion
  double epsilon_flavor(const std::vector<double> &k, unsigned int flavor) const;
  ///dispersion relation for a system with an additional phase
  double epsilon_phi(const std::vector<double> &k, const std::vector<double> &phi) const;

protected:


  std::vector<double> epsav() const { return epsav_; };
  std::vector<double> epssqav() const { return epssqav_; };
  double dispersion_moments(const unsigned int i) const { return dispersion_moments_[i]; }
  ///get epsilon(k) and the simpson integration weights on the BZ
  void compute_integration_weight_and_dispersion();
  ///initialize the vector phi with phases
  void setup_phase_vector(const alps::params &p);

  ///compute the center of gravity of the bands etc (used for high frequency analysis)
  void compute_dispersion_weights();

  ///reading input file
  template<typename GF> bool read_initial_G0(GF &G0_omega) const;

  ///compute G0 from an atomic solution
  virtual void compute_isolated_cluster_G0(cluster_nambu_matsubara_kspace_gf &G0_omega) const{ throw std::logic_error("this should be called in a subclass");}
  virtual void compute_isolated_cluster_G0(cluster_matsubara_kspace_gf &G0_omega) const{ throw std::logic_error("this should be called in a subclass");}
  ///compute G0 from a noninteracting solution
  virtual void compute_noninteracting_G0(cluster_nambu_matsubara_kspace_gf &G0_omega) const{ throw std::logic_error("this should be called in a subclass");};
  virtual void compute_noninteracting_G0(cluster_matsubara_kspace_gf &G0_omega) const{ throw std::logic_error("this should be called in a subclass");};

  ///chemical potential
  double mu_;
  double half_filling_;
  ///t is the basic energy unit of the hopping. t2 is the tprime. t_flavor is only used for mass imbalance
  double t_, t2_,t_flavor_;
  double beta_;
  int n_flavors_;
  int nfreq_;
  bool use_phase_;
  bool use_imbalance_;

  ///dispersion, appropriately binned for cluster point, lattice point (belonging to cluster), and spin (ignored in nambu case)
  boost::multi_array<double,3> epsilon_kckl_;
  std::vector<double> weight_kl_;
  std::vector<double> phi_;

  std::vector<double> epsav_;
  std::vector<double> epssqav_;
  std::vector<double> dispersion_moments_;
  std::string g0_input_path_;
  std::string initial_g0_;

};
class NormalStateClusterTransformer: public ClusterTransformer{
public:
  ///Constructor: a cluster transformer has no internal components, everything is delegated to the ClusterHelper
  NormalStateClusterTransformer(const alps::params& params): ClusterTransformer(params){}
  ///virtual destructor
  virtual ~NormalStateClusterTransformer() {}
  
  ///Use Dyson to get the self-energy
  void dyson_get_self_energy(cluster_matsubara_kspace_gf &sigma, const cluster_matsubara_kspace_gf &G, const cluster_matsubara_kspace_gf &G0) const;
  ///Use Dyson to get the G0
  void dyson_get_green0(cluster_matsubara_kspace_gf &G0, const cluster_matsubara_kspace_gf &sigma, const cluster_matsubara_kspace_gf &G) const;
  ///perform the k-sum to get G
  void matsubara_k_sum(cluster_matsubara_kspace_gf &G, const cluster_matsubara_kspace_gf &sigma) const;

  ///perform spin symmetrization of imaginary time k-space cluster Green's function
  void spin_symmetrize(cluster_itime_kspace_gf& G) const;
  ///perform spin symmetrization of matsubara k-space Green's function
  void spin_symmetrize(cluster_matsubara_kspace_gf& G) const;
  ///perform spin symmetrization of imaginary time k-space cluster Green's function
  void spin_symmetrize(cluster_itime_realspace_gf& G) const;
  ///perform spin symmetrization of matsubara k-space Green's function
  void spin_symmetrize(cluster_matsubara_realspace_gf& G) const;

  ///read G0 from file
  bool read_initial_G0 (cluster_matsubara_kspace_gf &G0_omega) const;
  ///compute G0 from an atomic solution
  void compute_isolated_cluster_G0(cluster_matsubara_kspace_gf &G0_omega) const;
  ///compute G0 from a noninteracting solution
  void compute_noninteracting_G0(cluster_matsubara_kspace_gf &G0_omega) const;

  void attach_bare_tail(cluster_matsubara_kspace_gf &G0) const;
  void attach_interacting_tail(cluster_matsubara_kspace_gf &G, double U, const std::vector<double> &densities) const;

  //Real space and k-space transformers for the Green's functions
  cluster_matsubara_realspace_gf transform_into_real_space(const cluster_matsubara_kspace_gf& G) const;
  cluster_matsubara_kspace_gf transform_into_k_space(const cluster_matsubara_realspace_gf& G) const;
  cluster_itime_realspace_gf transform_into_real_space(const cluster_itime_kspace_gf& G) const;
  cluster_itime_kspace_gf transform_into_k_space(const cluster_itime_realspace_gf& G) const;

  ///perform the DCA self-consistency cycle
  cluster_matsubara_kspace_gf operator()(const cluster_matsubara_kspace_gf &G_omega, const cluster_matsubara_kspace_gf &G0_omega) const;

private:
  //symmetrization in real momentum space - for the normal GF
  template<class GF,class TAIL> void symmetrize_helper_k_space(alps::gf::three_index_gf_with_tail<GF,TAIL> &G,
                                                        const std::pair<std::set<pair_class_type>::const_iterator,
                                                        std::set<pair_class_type>::const_iterator> pair_classes) const;
  template<class GF, class TAIL> void symmetrize_helper_real_space(alps::gf::four_index_gf_with_tail<GF,TAIL> &G,
                                                        const std::pair<std::set<pair_class_type>::const_iterator,
                                                        std::set<pair_class_type>::const_iterator> pair_classes) const;
public:
  template<class GF> void symmetrize_real_space(GF& G) const{return symmetrize_helper_real_space(G,pair_classes_real_space());}
  template<class GF> void symmetrize_k_space(GF& G) const{return symmetrize_helper_k_space(G,pair_classes_k_space());}
};

class DoubleCellClusterTransformer: public ClusterTransformer{
public:
  ///Constructor: a cluster transformer has no internal components, everything is delegated to the ClusterHelper
  DoubleCellClusterTransformer(const alps::params& params): ClusterTransformer(params),h_init_(params["dca.H_INIT"]), h_(params["dca.H"]){
    if(use_phase_) throw std::invalid_argument("use phase should not be used in the double unit cell framework.");
  }
  ///virtual destructor
  virtual ~DoubleCellClusterTransformer() {}

  ///read G0 from file
  bool read_initial_G0 (double_cell_matsubara_kspace_gf &G0_omega) const;

  ///compute G0 from an atomic solution
  void compute_isolated_cluster_G0(double_cell_matsubara_kspace_gf &G0_omega) const;
  ///compute G0 from a noninteracting solution
  void compute_noninteracting_G0(double_cell_matsubara_kspace_gf &G0_omega) const;

  ///the k-sum
  void matsubara_k_sum(double_cell_matsubara_kspace_gf &G, const double_cell_matsubara_kspace_gf &sigma, double h) const;
private:
  ///attach a tail to the noninteracting GF
  void attach_bare_tail(double_cell_matsubara_kspace_gf &G0, double h) const;

  ///initial magnetic field value
  const double h_init_;
  ///magnetic field
  const double h_;

};

class NambuClusterTransformer : public ClusterTransformer{
  ///Constructor: a cluster transformer has no internal components, everything is delegated to the ClusterHelper
public:
  NambuClusterTransformer(const alps::params& params);
  ///virtual destructor
  virtual ~NambuClusterTransformer(){};///perform the self-consistency transform. In need of a bit of refactoring; this is not nice.
  cluster_nambu_matsubara_kspace_gf operator()(const cluster_nambu_matsubara_kspace_gf & G_omega,
                                                      const cluster_nambu_matsubara_kspace_gf &G0_omega,
                                                      cluster_nambu_matsubara_kspace_gf &selfenergy,
                                                      double &mu, const double h, const double eta1, const double eta2,
                                                      const double beta, const bool use_selfenergy){ throw std::logic_error("a nambu k-space transform is not implemented for this transformer");}
  cluster_nambu_itime_kspace_gf symmetrize(const cluster_nambu_itime_kspace_gf& G) const{ throw std::logic_error("not implemented");}
  cluster_nambu_matsubara_kspace_gf symmetrize(const cluster_nambu_matsubara_kspace_gf& G) const{ throw std::logic_error("not implemented");}
  
  ///perform spin symmetrization of matsubara k-space Green's function
  void spin_symmetrize(cluster_nambu_matsubara_kspace_gf& G) const;

  ///read G0 from file
  bool read_initial_G0 (cluster_nambu_matsubara_kspace_gf &G0_omega) const;
  ///compute G0 from an atomic solution
  void compute_isolated_cluster_G0(cluster_nambu_matsubara_kspace_gf &G0_omega) const;
  ///compute G0 from a noninteracting solution
  void compute_noninteracting_G0(cluster_nambu_matsubara_kspace_gf &G0_omega) const;

  ///Use Dyson to get the self-energy
  void dyson_get_self_energy(cluster_nambu_matsubara_kspace_gf &sigma, const cluster_nambu_matsubara_kspace_gf &G, const cluster_nambu_matsubara_kspace_gf &G0) const;
  ///Use Dyson to get the G0
  void dyson_get_green0(cluster_nambu_matsubara_kspace_gf &G0, const cluster_nambu_matsubara_kspace_gf &sigma, const cluster_nambu_matsubara_kspace_gf &G) const;
  ///perform the k-sum to get G
  void matsubara_k_sum(cluster_nambu_matsubara_kspace_gf &G, const cluster_nambu_matsubara_kspace_gf &sigma) const;

  //Real space and k-space transformers for the Nambu Green's function
  virtual cluster_nambu_matsubara_realspace_gf transform_into_real_space(const cluster_nambu_matsubara_kspace_gf& G) const;
  virtual cluster_nambu_matsubara_kspace_gf transform_into_k_space(const cluster_nambu_matsubara_realspace_gf& G) const{ throw std::logic_error("not implemented");}
  virtual cluster_nambu_itime_realspace_gf transform_into_real_space(const cluster_nambu_itime_kspace_gf& G) const;
  virtual cluster_nambu_itime_kspace_gf transform_into_k_space(const cluster_nambu_itime_realspace_gf& G) const{ throw std::logic_error("not implemented");}

  void attach_bare_tail(cluster_nambu_matsubara_kspace_gf &G0) const;
  void attach_interacting_tail(cluster_nambu_matsubara_kspace_gf &G, double U, const std::vector<double> &densities) const;
  void symmetrize_dwave_kspace(cluster_nambu_matsubara_kspace_gf &G) const;
  void symmetrize_k_space(cluster_nambu_matsubara_kspace_gf &G) const{return symmetrize_dwave_kspace(G);}

private:
  ///superconducting bias field eta1, initial value
  double eta1_init_;
  ///superconducting bias field eta2, initial value
  double eta2_init_;
  ///superconducting bias field eta1
  double eta1_;
  ///superconducting bias field eta2
  double eta2_;
  ///buffered d-wave prefactor stored at the lattice momenta
  boost::multi_array<double,2>  dwave_prefactor_kckl_;
};

  

template<class GF,class TAIL> void NormalStateClusterTransformer::symmetrize_helper_k_space(alps::gf::three_index_gf_with_tail<GF,TAIL> &G,
                                                    const std::pair<std::set<pair_class_type>::const_iterator,
                                                    std::set<pair_class_type>::const_iterator> pair_classes) const
{
  std::set<pair_class_type>::const_iterator itc, itc_end;
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
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
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
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
private:
  std::complex<double> val[4];
};
