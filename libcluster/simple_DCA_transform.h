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
#include "cluster_transform.h"
#include "DCA_transform.h"

//does the DCA cluster transform for the paramagnetic phase (no AFM order)
class SimpleDCATransformer : public DCATransformer
{
public:
  SimpleDCATransformer(const alps::params& params) :
  DCATransformer(params),
  compute_filling_(params.defined("FILLING")),
  n_target_(compute_filling_?(double)(params["FILLING"]):-1)
  {
    
  }
  //self-consistency transform operator
  matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                        const matsubara_green_function_t &G0_omega,
                                        matsubara_green_function_t &selfenergy,
                                        double &mu, const double h, const double beta, const bool use_selfenergy);
  //legacy operator without self-energy
  matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                        const matsubara_green_function_t &G0_omega,
                                        double &mu, const double h, const double beta);
  //Hilbert transform operator for Nambu GF
  nambu_matsubara_green_function_t operator()(const nambu_matsubara_green_function_t &G_omega,
                                              const nambu_matsubara_green_function_t &G0_omega,
                                              nambu_matsubara_green_function_t &selfenergy,
                                              double &mu, const double h, const double eta1, const double eta2, const double beta, const bool use_selfenergy);
  
  //spin symmetrization functions
  itime_green_function_t symmetrize(const itime_green_function_t& G, const bool paramagnet) const
  { return spin_symmetrizer(G, paramagnet); }
  matsubara_green_function_t symmetrize(const matsubara_green_function_t& G, const bool paramagnet) const
  { return spin_symmetrizer(G, paramagnet); }
  matsubara_green_function_t transform_into_real_space(const matsubara_green_function_t& G) const
  {
    if(G.shape()!= diagonal) throw std::logic_error("transform of GF to real space: shape incorrect.");
    return transform_coordinates<std::complex<double>, Backward>(G);
  }
  matsubara_green_function_t transform_into_k_space(const matsubara_green_function_t& G) const
  {
    if(G.shape()!= nondiagonal) throw std::logic_error("transform of GF to kspace: shape incorrect.");
    return transform_coordinates<std::complex<double>, Forward>(G);
  }
  itime_green_function_t transform_into_real_space(const itime_green_function_t& G) const
  {
    if(G.shape()!= diagonal) throw std::logic_error("transform of GF to real space: shape incorrect.");
    return transform_coordinates<double, Backward>(G);
  }
  itime_green_function_t transform_into_k_space(const itime_green_function_t& G) const
  {
    if(G.shape()!= nondiagonal) throw std::logic_error("transform of GF to k space: shape incorrect.");
    return transform_coordinates<double, Forward>(G);
  }
  
  //spin symmetrization functions in Nambu space
  nambu_itime_green_function_t symmetrize(const nambu_itime_green_function_t& G, const bool paramagnet) const
  {
    std::cout<<"tried to symmetrize nambu itime."<<std::endl;
    return spin_symmetrizer(G, paramagnet);
  }
  nambu_matsubara_green_function_t symmetrize(const nambu_matsubara_green_function_t& G, const bool paramagnet) const
  {
    std::cout<<"tried to symmetrize nambu matsubara."<<std::endl;
    return spin_symmetrizer(G, paramagnet);
  }
  nambu_matsubara_green_function_t transform_into_real_space(const nambu_matsubara_green_function_t& G) const
  {
    return transform_coordinates<std::complex<double>, Backward>(G);
  }
  nambu_matsubara_green_function_t transform_into_k_space(const nambu_matsubara_green_function_t& G) const
  { return transform_coordinates<std::complex<double>, Forward>(G); }
  nambu_itime_green_function_t transform_into_real_space(const nambu_itime_green_function_t& G) const
  { return transform_coordinates<double, Backward>(G); }
  nambu_itime_green_function_t transform_into_k_space(const nambu_itime_green_function_t& G) const
  { return transform_coordinates<double, Forward>(G); }
  
  //initial G0
  virtual matsubara_green_function_t initial_G0      (const alps::params& parms, shape_t shape=diagonal);
  virtual nambu_matsubara_green_function_t initial_G0_nambu(const alps::params& parms);
  void construct_green_function(double mu, double beta, double eta1, double eta2, const std::vector<std::vector<std::vector<double> > > &epsilon_kckl, const std::vector<std::vector<double> > &dwave_prefactor_kckl, const std::vector<double> &weight_kl, const nambu_matsubara_green_function_t &selfenergy,nambu_matsubara_green_function_t &G_omega_lattice_kav);
  double compute_density(const nambu_matsubara_green_function_t &G_omega_lattice_kav, double beta);
  
private:
  template<typename T, class D>
  template <class T>
  green_function<T> spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const;
  
  //in nambu case: <double> very different from <complex<double> > -> no template
  nambu_green_function<std::complex<double> > spin_symmetrizer(const nambu_green_function<std::complex<double> >& G, const bool paramagnet) const;
  nambu_green_function<double> spin_symmetrizer(const nambu_green_function<double>& G, const bool paramagnet) const;
  
  matsubara_green_function_t fourpoint_green_function_inverse(const matsubara_green_function_t& G, const std::size_t nt) const;
  void compute_sigma_from_G0_and_G(
                                   const matsubara_green_function_t& G0_omega,
                                   const matsubara_green_function_t& G_omega,
                                   matsubara_green_function_t& selfenergy) const;
  void compute_integration_weight_and_dispersion(bool use_phase,std::vector<std::vector<std::vector<double> > > &epsilon_kckl, std::vector<double> &weight_kl, unsigned int nflavor) const;
  void compute_lattice_green_function(matsubara_green_function_t &G_omega_lattice, const matsubara_green_function_t &selfenergy, const double &mu, const double &beta, const double &h, std::vector<std::vector<std::vector<double> > > &epsilon_kckl, std::vector<double> &weight_kl) const;
  double compute_density(const matsubara_green_function_t &G_omega_lattice, const double &beta);
  void update_parameter_file(int current_it, double& mu);
  template<typename T>
  green_function<T> fourpoint_green_function_add(const green_function<T>& G1,
                                                 const green_function<T>& G2) const;
  template<typename T>
  green_function<T> fourpoint_green_function_subtract(const green_function<T>& G1,
                                                      const green_function<T>& G2) const;
  template<typename T>
  green_function<T> fourpoint_green_function_divide(const green_function<T>& G1,
                                                    const T d) const;
  
  struct Forward;
  struct Backward;
  
  ///should we tune the chemical potential to get the right density (compute_filling_=true) or should we work at constant mu?
  bool compute_filling_;
  ///target density. If compute_filling_ is true, we try to change mu to adjust to this.
  double n_target_;
};



