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

#include <set>

#include "boost/numeric/ublas/matrix.hpp"
#include "boost/tuple/tuple.hpp"

#include "alps/ngs/params.hpp"

#include "green_function.h"
#include "nambu_green_function.h"

#include "alpsgraph.h"
#include "clusterhelper.h"


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


//The DCA transormer does the Fourier transforms, dispersions, etc that are needed for the DCA.
class DCATransformer : public ClusterTransformer
{
public:
  DCATransformer(const alps::params& params);
  virtual ~DCATransformer(){} //pro forma virtual destructor
  
  //symmetrization in real momentum space - for the normal GF
  template<typename T> green_function<T> symmetrize_real_space(const green_function<T>& G) const
  { return symmetrize_helper(G, pair_classes_real_space()); }
  template<typename T> green_function<T> symmetrize_k_space(const green_function<T>& G) const
  { return symmetrize_helper(G, pair_classes_k_space()); }
  
  //symmetrization in real and momentum space - for the NAMBU GF
  template<typename T> nambu_green_function<T> symmetrize_real_space(const nambu_green_function<T>& G) const
  { return symmetrize_helper(G, pair_classes_real_space()); }
  template<typename T> nambu_green_function<T> symmetrize_k_space(const nambu_green_function<T>& G) const
  { return symmetrize_helper(G, pair_classes_k_space()); }
  const Matrix<int>::Type &symmetry_table_k_space() const{return ClusterHelper::symmetry_table_k_space();}
  const std::set<pair_class_type> &pair_class_set_k_space() const{return ClusterHelper::pair_class_set_k_space();}
  
  
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
  
  virtual matsubara_green_function_t initial_G0      (const alps::params& parms, shape_t shape=diagonal);
  //Hilbert transform operator
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                                const matsubara_green_function_t &G0_omega,
                                                double &mu, const double h, const double beta)=0;
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                                const matsubara_green_function_t &G0_omega,
                                                matsubara_green_function_t &selfenergy,
                                                double &mu, const double h, const double beta, const bool use_selfenergy)=0;
  //Hilbert transform operator for Nambu GF
  virtual nambu_matsubara_green_function_t operator()(const nambu_matsubara_green_function_t &G_omega,
                                                      const nambu_matsubara_green_function_t &G0_omega,
                                                      nambu_matsubara_green_function_t &selfenergy,
                                                      double &mu, const double h, const double eta1, const double eta2, const double beta, const bool use_selfenergy)=0;
  
  int* momentum_translation() const;
  std::vector<double> epsav() const { return epsav_; };
  std::vector<double> epssqav() const { return epssqav_; };
  double dispersion_moments(const unsigned int i) const { return dispersion_moments_[i]; }
  double bandwidth() const { return bandwidth_; };
  std::vector<double> tmatrix() const;
  virtual void write_hdf5(alps::hdf5::archive &ar, const std::string &p) const;
  
  //Generate a header string that can be printed at the beginning of a Greens function file to label the columns
  std::string generate_header_kspace(const matsubara_green_function_t & G) const;
  std::string generate_header_realspace(const matsubara_green_function_t & G) const;
  std::string generate_header_kspace(const itime_green_function_t & G) const;
  std::string generate_header_realspace(const itime_green_function_t & G) const;
  void print_momenta(const char *filename, const matsubara_green_function_t & G) const;
  
  double epsilon(const std::vector<double> &k) const
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
    if(t3_!=0.) {
      switch(dimension()) {
        case 1 :
          eps -= 2*t3_*cos(k[0]*3.);
          break;
        case 2 :
          break;
        case 3 :
          eps -= 2*t3_*cos(k[0]+k[1]+k[2]);
          eps -= 2*t3_*cos(k[0]-k[1]+k[2]);
          eps -= 2*t3_*cos(k[0]+k[1]-k[2]);
          eps -= 2*t3_*cos(k[0]-k[1]-k[2]);
          break;
      }
    }
    if (t4_!=0.) {
      switch(dimension()) {
        case 1 :
          eps -= 2*t4_*cos(k[0]*4.);
          break;
        case 2 :
          break;
        case 3 :
          eps -= 2*t4_*cos(k[0]*2.);
          eps -= 2*t4_*cos(k[1]*2.);
          eps -= 2*t4_*cos(k[2]*2.);
          break;
      }
    }
    return eps*bandwidth_scale_factor_;
  }
  //dispersion, evaluated with a phase shift of phi
  double epsilon_phi(const std::vector<double> &k, const std::vector<double> &phi) const;
  std::vector<double> phase_vector() const{
    std::vector<double> phi(dimension());
    for(int i=0;i<dimension();++i){
      std::stringstream phase_name; phase_name<<"PHASE_"<<i;
      phi[i]=params_[phase_name.str()];
      phi[i] = phi[i]/periodicity(i);
    }
    return phi;
  }
  
  //dispersion, evaluated with flavor dependent t_
  double epsilon_flavor(const std::vector<double> &k, unsigned int flavor) const
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
    
    return eps*bandwidth_scale_factor_;
  }
  
  
private:
  template<typename T>
  green_function<T> symmetrize_helper(const green_function<T>& G,
                                      const std::pair<std::set<pair_class_type>::const_iterator,
                                      std::set<pair_class_type>::const_iterator> pair_classes) const;
  template<typename T>
  nambu_green_function<T> symmetrize_helper(const nambu_green_function<T>& G,
                                            const std::pair<std::set<pair_class_type>::const_iterator,
                                            std::set<pair_class_type>::const_iterator> pair_classes) const;
  const double t_, t2_, t3_, t4_;
  const double t_flavor_;
  double bandwidth_;
  double bandwidth_scale_factor_;
  double half_filling_;
  std::vector<double> epsav_;
  std::vector<double> epssqav_;
  std::vector<double> dispersion_moments_;
  
protected:
  const alps::params params_;
};


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
  green_function<T> transform_coordinates(const       green_function<T>& G) const;
  template<typename T, class D>
  nambu_green_function<T> transform_coordinates(const nambu_green_function<T>& G) const;
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



class AntiFerroDCATransformer : public DCATransformer
{
public:
  AntiFerroDCATransformer(const alps::params& params) : DCATransformer(params) {}
  
  //self-consistency transform operator
  matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                        const matsubara_green_function_t &G0_omega,
                                        matsubara_green_function_t &selfenergy,
                                        double &mu, const double h, const double beta, const bool use_selfenergy);
  matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                        const matsubara_green_function_t &G0_omega,
                                        double &mu, const double h, const double beta);
  nambu_matsubara_green_function_t operator()(const nambu_matsubara_green_function_t &G_omega,
                                              const nambu_matsubara_green_function_t &G0_omega,
                                              nambu_matsubara_green_function_t &selfenergy,
                                              double &mu, const double h, const double eta1, const double eta2, const double beta, const bool use_selfenergy);
  itime_green_function_t symmetrize(const itime_green_function_t& G, const bool paramagnet) const
  { throw std::runtime_error("spin symmetrization in AFM not implemented for itime.");}//return spin_symmetrizer(G, paramagnet); }
  matsubara_green_function_t symmetrize(const matsubara_green_function_t& G, const bool paramagnet) const
  { return spin_symmetrizer(G, paramagnet); }
  matsubara_green_function_t transform_into_real_space(const matsubara_green_function_t& G) const
  {
    if(G.shape()!= blockdiagonal) throw std::logic_error("transform of GF to real space: shape incorrect.");
    return transform_coordinates<std::complex<double>, Backward>(G);
  }
  matsubara_green_function_t transform_into_k_space(const matsubara_green_function_t& G) const
  {
    if(G.shape()!= nondiagonal) throw std::logic_error("transform of GF to k space: shape incorrect.");
    return transform_coordinates<std::complex<double>, Forward>(G);
  }
  itime_green_function_t transform_into_real_space(const itime_green_function_t& G) const
  {
    if(G.shape()!= blockdiagonal) throw std::logic_error("transform of GF to real space: shape incorrect.");
    return transform_coordinates<double, Backward>(G);
  }
  itime_green_function_t transform_into_k_space(const itime_green_function_t& G) const
  {
    if(G.shape()!= nondiagonal) throw std::logic_error("transform of GF to k space: shape incorrect.");
    return transform_coordinates<double, Forward>(G);
  }
  
  nambu_matsubara_green_function_t transform_into_real_space(const nambu_matsubara_green_function_t& G) const {throw std::logic_error("AFM NAMBU not implemented.");}
  //{ return transform_coordinates<std::complex<double>, Backward>(G); }
  nambu_matsubara_green_function_t transform_into_k_space(const nambu_matsubara_green_function_t& G) const {throw std::logic_error("AFM NAMBU not implemented.");}
  //{ return transform_coordinates<std::complex<double>, Forward>(G); }
  nambu_itime_green_function_t transform_into_real_space(const nambu_itime_green_function_t& G) const {throw std::logic_error("AFM NAMBU not implemented.");}
  //{ return transform_coordinates<double, Backward>(G); }
  nambu_itime_green_function_t transform_into_k_space(const nambu_itime_green_function_t& G) const {throw std::logic_error("AFM NAMBU not implemented.");}
  //{ return transform_coordinates<double, Forward>(G); }
  
  
  virtual matsubara_green_function_t initial_G0(const alps::params& parms, shape_t shape=blockdiagonal);
  static matsubara_green_function_t green_function_inverse(const matsubara_green_function_t& G);
  static matsubara_green_function_t green_function_multiply(const matsubara_green_function_t& G1, const matsubara_green_function_t& G2);
  static matsubara_green_function_t green_function_add(const matsubara_green_function_t& G1, const matsubara_green_function_t& G2);
  static matsubara_green_function_t green_function_subtract(const matsubara_green_function_t& G1, const matsubara_green_function_t& G2);
  static matsubara_green_function_t green_function_divide(const matsubara_green_function_t& G1, const std::complex<double> d);
  
  std::map<double, ClusterHelper::Matrix<std::complex<double> >::Type>
  local_green_function_from_selfenergies(const std::map<double,
                                         std::vector<ClusterHelper::Matrix<std::complex<double> >::Type> >& selfenergies,
                                         const double mu, const double delta);
  
private:
  matsubara_green_function_t spin_symmetrizer(const matsubara_green_function_t& G, const bool paramagnet) const;
  template<typename T, class D>
  green_function<T> transform_coordinates(const green_function<T>& G) const;
  
  struct Forward;
  struct Backward;
  
};

