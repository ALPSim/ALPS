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


//The Cluster transormer presents a framework for operator(), which does the Hilbert transform..
class ClusterTransformer : protected ClusterHelper
{
public:
  ClusterTransformer(const alps::params& params) :
  ClusterHelper(params) {
    //std::cerr<<"done building cluster transformer"<<std::endl;
  }
  static ClusterTransformer* generate_transformer(const alps::params& parms);
  virtual ~ClusterTransformer() {} //pro forma virtual destructor.
  
  ///perform a Hilber transform: given a G and a G0, produce a G0 for the next iteration.
  virtual       matsubara_green_function_t operator()(const matsubara_green_function_t & G_omega,
                                                      const matsubara_green_function_t &G0_omega,
                                                      double &mu, const double h,
                                                      const double beta)=0;
  virtual nambu_matsubara_green_function_t operator()(const nambu_matsubara_green_function_t & G_omega,
                                                      const nambu_matsubara_green_function_t &G0_omega,
                                                      nambu_matsubara_green_function_t &selfenergy,
                                                      double &mu, const double h, const double eta1, const double eta2,
                                                      const double beta, const bool use_selfenergy)=0;
  ///perform symmetrizaton in imaginary time
  virtual       itime_green_function_t symmetrize(const       itime_green_function_t& G, const bool paramagnet) const=0;
  virtual nambu_itime_green_function_t symmetrize(const nambu_itime_green_function_t& G, const bool paramagnet) const{ throw std::logic_error("not implemented");}
  ///perform symmetrizaton in matsubara frequencies time
  virtual       matsubara_green_function_t symmetrize(const       matsubara_green_function_t& G, const bool paramagnet) const=0;
  virtual nambu_matsubara_green_function_t symmetrize(const nambu_matsubara_green_function_t& G, const bool paramagnet) const{ throw std::logic_error("not implemented");};
  ///compute an initial G0.
  virtual       matsubara_green_function_t initial_G0      (const alps::params& parms, shape_t shape=nondiagonal);
  
  virtual nambu_matsubara_green_function_t initial_G0_nambu(const alps::params& parms);
  virtual std::vector<double> epsav() const=0;
  virtual std::vector<double> epssqav() const=0;
  nambu_matsubara_green_function_t initial_G_nambu(const alps::params& parms);
  
  //Real space and k-space transformers for the Green's functions
  virtual matsubara_green_function_t transform_into_real_space(const matsubara_green_function_t& G) const=0;
  virtual matsubara_green_function_t transform_into_k_space(const matsubara_green_function_t& G) const=0;
  virtual itime_green_function_t transform_into_real_space(const itime_green_function_t& G) const=0;
  virtual itime_green_function_t transform_into_k_space(const itime_green_function_t& G) const=0;
  //Real space and k-space transformers for the Nambu Green's function
  virtual nambu_matsubara_green_function_t transform_into_real_space(const nambu_matsubara_green_function_t& G) const=0;
  virtual nambu_matsubara_green_function_t transform_into_k_space(const nambu_matsubara_green_function_t& G) const=0;
  virtual nambu_itime_green_function_t transform_into_real_space(const nambu_itime_green_function_t& G) const=0;
  virtual nambu_itime_green_function_t transform_into_k_space(const nambu_itime_green_function_t& G) const=0;
  
  std::complex<double>* exp_iksite(const int sign) const;
  
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
};



