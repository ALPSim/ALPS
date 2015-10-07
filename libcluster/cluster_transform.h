/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2007      by Sebastian Fuchs <fuchs@comp-phys.org>
 * Copyright (C) 2007-1013 by Emanuel Gull <egull@umich.edu>
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

//The Cluster transformer presents a framework for operator(), which does the Hilbert transform..
class NormalStateClusterTransformer : public ClusterHelper
{
public:
  ///Constructor: a cluster transformer has no internal components, everything is delegated to the ClusterHelper
  NormalStateClusterTransformer(const alps::params& params) : ClusterHelper(params) {  }
  ///virtual destructor
  virtual ~NormalStateClusterTransformer() {}
  
  ///perform the self-consistency transform: given a G and a G0, produce a G0 for the next iteration.
  virtual       cluster_matsubara_kspace_gf operator()(const cluster_matsubara_kspace_gf & G_omega,
                                                      const cluster_matsubara_kspace_gf &G0_omega,
                                                      double &mu, const double h,
                                                      const double beta){ throw std::logic_error("a normal-state k-space transform is not implemented for this transformer");}
  ///perform spin symmetrization of imaginary time k-space cluster Green's function
  void spin_symmetrize(cluster_itime_kspace_gf& G) const;
  ///perform spin symmetrization of matsubara k-space Green's function
  void spin_symmetrize(cluster_matsubara_kspace_gf& G) const;
  ///perform spin symmetrization of imaginary time k-space cluster Green's function
  void spin_symmetrize(cluster_itime_realspace_gf& G) const;
  ///perform spin symmetrization of matsubara k-space Green's function
  void spin_symmetrize(cluster_matsubara_realspace_gf& G) const;

  //symmetrization in real momentum space - for the normal GF
  template<class GF,class TAIL> void symmetrize_helper_k_space(alps::gf::three_index_gf_with_tail<GF,TAIL> &G,
                                                        const std::pair<std::set<pair_class_type>::const_iterator,
                                                        std::set<pair_class_type>::const_iterator> pair_classes) const;
  template<class GF, class TAIL> void symmetrize_helper_real_space(alps::gf::four_index_gf_with_tail<GF,TAIL> &G,
                                                        const std::pair<std::set<pair_class_type>::const_iterator,
                                                        std::set<pair_class_type>::const_iterator> pair_classes) const;
  template<class GF> void symmetrize_real_space(GF& G) const{return symmetrize_helper_real_space(G,pair_classes_real_space());}
  template<class GF> void symmetrize_k_space(GF& G) const{return symmetrize_helper_k_space(G,pair_classes_k_space());};


  ///compute an initial G0.
  virtual cluster_matsubara_kspace_gf initial_G0 (const alps::params& parms);

  //->probably promote to super class
  virtual std::vector<double> epsav() const{ throw std::logic_error("not implemented");};
  virtual std::vector<double> epssqav() const{ throw std::logic_error("not implemented");};

  //Real space and k-space transformers for the Green's functions
  cluster_matsubara_realspace_gf transform_into_real_space(const cluster_matsubara_kspace_gf& G) const;
  cluster_matsubara_kspace_gf transform_into_k_space(const cluster_matsubara_realspace_gf& G) const;
  cluster_itime_realspace_gf transform_into_real_space(const cluster_itime_kspace_gf& G) const;
  cluster_itime_kspace_gf transform_into_k_space(const cluster_itime_realspace_gf& G) const;

  alps::gf::momentum_index_mesh momentum_mesh() const;
  alps::gf::real_space_index_mesh real_space_mesh() const;

};

class NambuClusterTransformer : protected ClusterHelper{
  ///Constructor: a cluster transformer has no internal components, everything is delegated to the ClusterHelper
  NambuClusterTransformer(const alps::params& params) : ClusterHelper(params) {  }
  ///virtual destructor
  virtual ~NambuClusterTransformer() {}///perform the self-consistency transform. In need of a bit of refactoring; this is not nice.
  virtual cluster_nambu_matsubara_kspace_gf operator()(const cluster_nambu_matsubara_kspace_gf & G_omega,
                                                      const cluster_nambu_matsubara_kspace_gf &G0_omega,
                                                      cluster_nambu_matsubara_kspace_gf &selfenergy,
                                                      double &mu, const double h, const double eta1, const double eta2,
                                                      const double beta, const bool use_selfenergy){ throw std::logic_error("a nambu k-space transform is not implemented for this transformer");}
  virtual cluster_nambu_itime_kspace_gf symmetrize(const cluster_nambu_itime_kspace_gf& G) const{ throw std::logic_error("not implemented");}
  virtual cluster_nambu_matsubara_kspace_gf symmetrize(const cluster_nambu_matsubara_kspace_gf& G) const{ throw std::logic_error("not implemented");}
  
  virtual cluster_nambu_matsubara_kspace_gf initial_G0_nambu(const alps::params& parms);

  cluster_nambu_matsubara_kspace_gf initial_G_nambu(const alps::params& parms);

  //Real space and k-space transformers for the Nambu Green's function
  virtual cluster_nambu_matsubara_realspace_gf transform_into_real_space(const cluster_nambu_matsubara_kspace_gf& G) const{ throw std::logic_error("not implemented");};
  virtual cluster_nambu_matsubara_kspace_gf transform_into_k_space(const cluster_nambu_matsubara_realspace_gf& G) const{ throw std::logic_error("not implemented");};
  virtual cluster_nambu_itime_realspace_gf transform_into_real_space(const cluster_nambu_itime_kspace_gf& G) const{ throw std::logic_error("not implemented");};
  virtual cluster_nambu_itime_kspace_gf transform_into_k_space(const cluster_nambu_itime_realspace_gf& G) const{ throw std::logic_error("not implemented");};
//  template<typename T, class D>
//  nambu_green_function<T> transform_coordinates(const nambu_green_function<T>& G) const;

};

  
 /* std::complex<double>* exp_iksite(const int sign) const;
  
  template<typename T>
  T interpolate(std::vector<double> k, const std::vector<T>& values, const bool antiferromagnet) const {
    return ClusterHelper::interpolate(k, values, antiferromagnet);
  }
  
  unsigned int dimension() const { return ClusterHelper::dimension(); }
  double cluster_momentum(const unsigned int i, const unsigned int j) const { return ClusterHelper::cluster_momentum(i,j); }
  double lattice_momentum(const unsigned int i, const unsigned int j) const { return ClusterHelper::lattice_momentum(i,j); }
  double reciprocal_basis_vector(const unsigned int i, const unsigned int j) const { return ClusterHelper::reciprocal_basis_vector(i,j); }
  int periodicity(const unsigned int dim) const { return ClusterHelper::periodicity(dim); }
  unsigned int n_lattice_momenta() const { return ClusterHelper::n_lattice_momenta(); }
  double site_distance(const unsigned int i, const unsigned int j) const { return ClusterHelper::site_distance(i, j); }
  virtual void write_hdf5(alps::hdf5::archive &ar, const std::string &p) const;
};*/


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


