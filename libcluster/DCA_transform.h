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




//The DCA transormer does the Fourier transforms, dispersions, etc that are needed for the DCA.
class DCATransformer : public ClusterTransformer
{
public:
  typedef std::pair<int, int> pair_type;
  typedef std::set<pair_type> pair_class_type;
  typedef std::set<pair_class_type> pair_class_set_type;

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
  const pair_class_set_type &pair_class_set_k_space() const{return ClusterHelper::pair_class_set_k_space();}
  
  
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


