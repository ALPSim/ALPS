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

#ifndef ALPS_CLUSTER_H
#define ALPS_CLUSTER_H

#include "boost/numeric/bindings/detail/config/fortran.hpp"
#include "boost/numeric/bindings/ublas.hpp"
#include "boost/numeric/bindings/lapack/driver/gesv.hpp"

#include "alps/lattice/graph_helper.h"
//#include "alps/alea.h"

#include "types.h"
#include "green_function.h"
#include "nambu_green_function.h"
#include "hilberttransformer.h"

//forward declarations
namespace alps{
  template <class T> class SimpleObservableEvaluator;
  typedef class SimpleObservableEvaluator<double> RealObsevaluator;
}
typedef boost::numeric::ublas::matrix<std::complex<double>,boost::numeric::ublas::column_major> complex_dense_matrix;


class ClusterHelper: private alps::graph_helper<alps::graph_helper<>::graph_type>
{
public:
  
  template <typename T> struct Matrix
  {
    typedef boost::numeric::ublas::matrix<T, boost::numeric::ublas::column_major> Type;
    
    static Type inverse(Type M);
    
    static Type multiply(Type M1, Type M2) {
      return boost::numeric::ublas::prec_prod(M1, M2);
    }
  };
  
  typedef alps::graph_helper<alps::graph_helper<>::graph_type>::vector_type vector_type;
  typedef std::set<site_t> site_class_type;
  typedef std::pair<site_t, site_t> pair_type;
  typedef std::set<pair_type> pair_class_type;
  typedef std::set<site_class_type> site_class_set_type;
  typedef std::set<pair_class_type> pair_class_set_type;
  typedef std::pair<site_class_set_type::const_iterator, site_class_set_type::const_iterator> site_class_set_itpair_type;
  typedef std::pair<pair_class_set_type::const_iterator, pair_class_set_type::const_iterator> pair_class_set_itpair_type;
  
  ClusterHelper(const alps::Parameters& params);
  
  double cluster_coordinate(const unsigned int i, const unsigned int j) const { return cluster_coordinates_(i,j); }
  double cluster_momentum(const unsigned int i, const unsigned int j) const { return cluster_momenta_(i,j); }
  double lattice_momentum(const unsigned int i, const unsigned int j) const { return lattice_momenta_(i,j); }
  bool bipartite() const { return bipartite_; }
  unsigned int n_cluster_sites() const { return n_cluster_sites_; }
  unsigned int n_cluster_cells() const { return alps::graph_helper<alps::graph_helper<>::graph_type>::volume(); }
  unsigned int n_lattice_momenta() const { return lattice_momenta_.size1(); }
  unsigned int dimension() const { return alps::graph_helper<alps::graph_helper<>::graph_type>::dimension(); }
  unsigned int cluster_momentum_approximation(vector_type k) const;
  
  template<typename T>
  T interpolate(vector_type k, const std::vector<T>& values, const bool antiferromagnet) const;
  double site_distance(const unsigned int i, const unsigned int j) const;
  
  site_class_set_itpair_type site_classes_real_space() const
  { return site_class_set_itpair_type (site_class_set_real_space_.begin(), site_class_set_real_space_.end()); }
  pair_class_set_itpair_type pair_classes_real_space() const
  { return pair_class_set_itpair_type(pair_class_set_real_space_.begin(), pair_class_set_real_space_.end()); }
  site_class_set_itpair_type site_classes_k_space() const
  { return site_class_set_itpair_type(site_class_set_k_space_.begin(), site_class_set_k_space_.end()); }
  pair_class_set_itpair_type pair_classes_k_space() const
  { return pair_class_set_itpair_type(pair_class_set_k_space_.begin(), pair_class_set_k_space_.end()); }
  const Matrix<site_t>::Type &symmetry_table_k_space() const{return symmetry_table_k_space_;}
  
private:
  void setup_tables(const alps::Parameters& params);
  void symmetry_analysis(const Matrix<double>::Type& basis,
                         const Matrix<double>::Type& basis_inv);
  void ladder_symmetry_analysis();
  std::pair<site_class_set_type, pair_class_set_type> setup_classes(const Matrix<site_t>::Type& symmetry_table) const;
  bool inside_cluster(vector_type v, const Matrix<double>::Type& b) const;
  bool inside_reduced_brillouin_zone(vector_type v) const;
  bool append_column_no_duplicates(Matrix<double>::Type& m, const vector_type& v) const;
  bool append_column_no_duplicates(Matrix<site_t>::Type& m, std::map<site_t, site_t>& smap) const;
  vector_type map_back_into_brillouin_zone(vector_type v) const;
  vector_type map_back_into_cluster(vector_type v,
                                    const Matrix<double>::Type& B_inv, const Matrix<double>::Type& B) const;
  unsigned int identify_vertex(const vector_type& v, const Matrix<double>::Type& coordinate_matrix) const;
  unsigned int identify_sublattice(const Matrix<double>::Type& coordinate_matrix, const unsigned int i) const;
  double round(const double d) const;
  vector_type round(vector_type v) const;
  bool almost_equal(const double d1, const double d2) const { return fabs(d1-d2)<epsilon_; }
  vector_type matrix_vector_mult(const Matrix<double>::Type& M, const vector_type& v) const;
  bool write(const alps::Parameters &parms) const;
  bool read(const alps::Parameters &parms);
  void write_site_class_set(alps::hdf5::archive &ar, const site_class_set_type &set, const std::string &name) const;
  void write_pair_class_set(alps::hdf5::archive &ar, const pair_class_set_type &set, const std::string &name) const;
  void read_site_class_set(alps::hdf5::archive &ar, site_class_set_type &set, const std::string &name);
  void read_pair_class_set(alps::hdf5::archive &ar, pair_class_set_type &set, const std::string &name);
  
  
  Matrix<double>::Type cluster_coordinates_;
  Matrix<double>::Type cluster_momenta_;
  Matrix<double>::Type lattice_momenta_;
  Matrix<double>::Type reciprocal_basis_;
  Matrix<double>::Type reciprocal_basis_inv_;
  Matrix<double>::Type rbzone_boundaries_;
  site_class_set_type site_class_set_real_space_;
  pair_class_set_type pair_class_set_real_space_;
  site_class_set_type site_class_set_k_space_;
  pair_class_set_type pair_class_set_k_space_;
  Matrix<site_t>::Type symmetry_table_k_space_;
  Matrix<site_t>::Type symmetry_table_real_space_;
  
  unsigned int n_cluster_sites_;
  bool bipartite_;
  const double epsilon_;
  std::string lattice_momenta_file_name_;
};


//The Cluster transormer presents a framework for operator(), which does the Hilbert transform..
class ClusterTransformer : protected ClusterHelper
{
public:
  ClusterTransformer(const alps::Parameters& params) :
  ClusterHelper(params) {}
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
  virtual       matsubara_green_function_t initial_G0      (const alps::Parameters& parms, shape_t shape=nondiagonal);
  
  virtual nambu_matsubara_green_function_t initial_G0_nambu(const alps::Parameters& parms);
  nambu_matsubara_green_function_t initial_G_nambu(const alps::Parameters& parms);
  
  template<typename T>
  T interpolate(vector_type k, const std::vector<T>& values, const bool antiferromagnet) const {
    return ClusterHelper::interpolate(k, values, antiferromagnet);
  }
  
  unsigned int dimension() const { return ClusterHelper::dimension(); }
  double site_distance(const unsigned int i, const unsigned int j) const { return ClusterHelper::site_distance(i, j); }
};


//The DCA transormer does the Fourier transforms, dispersions, etc that are needed for the DCA.
class DCATransformer : public ClusterTransformer
{
public:
  DCATransformer(const alps::Parameters& params);
  virtual ~DCATransformer(){} //pro forma virtual destructor
  static void generate_transformer(const alps::Parameters& parms, boost::shared_ptr<DCATransformer> &trans_ptr);
  
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
  const Matrix<site_t>::Type &symmetry_table_k_space() const{return ClusterHelper::symmetry_table_k_space();}
  
  
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
  
  virtual matsubara_green_function_t transform_fourpoint(const matsubara_green_function_t& G_omega_fourpoint,
                                                         const matsubara_green_function_t& G_omega,
                                                         const matsubara_green_function_t& G0_omega,
                                                         const double mu, const double h, const double beta,
                                                         const unsigned int n_matsubara) const =0;
  virtual alps::RealObsevaluator kinetic_energy(const std::vector<std::vector<std::vector<alps::RealObsevaluator> > >& selfenergies,
                                                const double beta, const double mu,
                                                const bool error_propagation, const bool paramagnet)=0;
  virtual matsubara_green_function_t initial_G0      (const alps::Parameters& parms, shape_t shape=diagonal);
  //Hilbert transform operator
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                                const matsubara_green_function_t &G0_omega,
                                                double &mu, const double h, const double beta)=0;
  //Hilbert transform operator for Nambu GF
  virtual nambu_matsubara_green_function_t operator()(const nambu_matsubara_green_function_t &G_omega,
                                                      const nambu_matsubara_green_function_t &G0_omega,
                                                      nambu_matsubara_green_function_t &selfenergy,
                                                      double &mu, const double h, const double eta1, const double eta2, const double beta, const bool use_selfenergy)=0;
  
  std::complex<double>* exp_iksite(const int sign) const;
  int* momentum_translation() const;
  std::vector<double> epsav() const { return epsav_; };
  std::vector<double> epssqav() const { return epssqav_; };
  double dispersion_moments(const unsigned int i) const { return dispersion_moments_[i]; }
  double bandwidth() const { return bandwidth_; };
  
  double epsilon(const vector_type k) const
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
  double cluster_momentum(const unsigned int i, const unsigned int j) const { return ClusterHelper::cluster_momentum(i,j); }
  
  
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
  double bandwidth_;
  double bandwidth_scale_factor_;
  double half_filling_;
  std::vector<double> epsav_;
  std::vector<double> epssqav_;
  std::vector<double> dispersion_moments_;
};


//does the DCA cluster transform for the paramagnetic phase (no AFM order)
class SimpleDCATransformer : public DCATransformer
{
public:
  SimpleDCATransformer(const alps::Parameters& params) :
  DCATransformer(params),
  compute_filling_(params.defined("FILLING")),
  n_target_(compute_filling_?(double)(params["FILLING"]):-1)
  {
    params_ptr_=&params;
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
  { return transform_coordinates<std::complex<double>, Backward>(G); }
  matsubara_green_function_t transform_into_k_space(const matsubara_green_function_t& G) const
  { return transform_coordinates<std::complex<double>, Forward>(G); }
  itime_green_function_t transform_into_real_space(const itime_green_function_t& G) const
  { return transform_coordinates<double, Backward>(G); }
  itime_green_function_t transform_into_k_space(const itime_green_function_t& G) const
  { return transform_coordinates<double, Forward>(G); }
  
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
  virtual matsubara_green_function_t initial_G0      (const alps::Parameters& parms, shape_t shape=diagonal);
  virtual nambu_matsubara_green_function_t initial_G0_nambu(const alps::Parameters& parms);
  matsubara_green_function_t transform_fourpoint(const matsubara_green_function_t& G_omega_foupoint,
                                                 const matsubara_green_function_t &G_omega,
                                                 const matsubara_green_function_t &G0_omega,
                                                 const double mu, const double h, const double beta,
                                                 const unsigned int n_matsubara) const;
  alps::RealObsevaluator kinetic_energy(const std::vector<std::vector<std::vector<alps::RealObsevaluator> > >& selfenergies,
                                        const double beta, const double mu,
                                        const bool error_propagation, const bool paramagnet);
  std::map<double, std::complex<double> > local_green_function_from_selfenergies(const std::map<double,
                                                                                 std::vector<std::complex<double> > >& selfenergies,
                                                                                 const double mu, const double delta);
  void construct_green_function(double mu, double beta, double eta1, double eta2, const std::vector<std::vector<double> > &epsilon_kckl, const std::vector<std::vector<double> > &dwave_prefactor_kckl, const std::vector<double> &weight_kl, const nambu_matsubara_green_function_t &selfenergy,nambu_matsubara_green_function_t &G_omega_lattice_kav);
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
  
  matsubara_green_function_t fourpoint_green_function_inverse(const matsubara_green_function_t& G, const unsigned int nt) const;
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
  
  bool compute_filling_;
  double n_target_;
  const alps::Parameters *params_ptr_;
};



class AntiFerroDCATransformer : public DCATransformer
{
public:
  AntiFerroDCATransformer(const alps::Parameters& params) : DCATransformer(params) {}
  
  matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                        const matsubara_green_function_t &G0_omega,
                                        double &mu, const double h, const double beta);
  nambu_matsubara_green_function_t operator()(const nambu_matsubara_green_function_t &G_omega,
                                              const nambu_matsubara_green_function_t &G0_omega,
                                              nambu_matsubara_green_function_t &selfenergy,
                                              double &mu, const double h, const double eta1, const double eta2, const double beta, const bool use_selfenergy);
  itime_green_function_t symmetrize(const itime_green_function_t& G, const bool paramagnet) const
  { return spin_symmetrizer(G, paramagnet); }
  matsubara_green_function_t symmetrize(const matsubara_green_function_t& G, const bool paramagnet) const
  { return spin_symmetrizer(G, paramagnet); }
  matsubara_green_function_t transform_into_real_space(const matsubara_green_function_t& G) const
  { return transform_coordinates<std::complex<double>, Backward>(G); }
  matsubara_green_function_t transform_into_k_space(const matsubara_green_function_t& G) const
  { return transform_coordinates<std::complex<double>, Forward>(G); }
  itime_green_function_t transform_into_real_space(const itime_green_function_t& G) const
  { return transform_coordinates<double, Backward>(G); }
  itime_green_function_t transform_into_k_space(const itime_green_function_t& G) const
  { return transform_coordinates<double, Forward>(G); }
  
  nambu_matsubara_green_function_t transform_into_real_space(const nambu_matsubara_green_function_t& G) const {throw std::logic_error("AFM NAMBU not implemented.");}
  //{ return transform_coordinates<std::complex<double>, Backward>(G); }
  nambu_matsubara_green_function_t transform_into_k_space(const nambu_matsubara_green_function_t& G) const {throw std::logic_error("AFM NAMBU not implemented.");}
  //{ return transform_coordinates<std::complex<double>, Forward>(G); }
  nambu_itime_green_function_t transform_into_real_space(const nambu_itime_green_function_t& G) const {throw std::logic_error("AFM NAMBU not implemented.");}
  //{ return transform_coordinates<double, Backward>(G); }
  nambu_itime_green_function_t transform_into_k_space(const nambu_itime_green_function_t& G) const {throw std::logic_error("AFM NAMBU not implemented.");}
  //{ return transform_coordinates<double, Forward>(G); }
  
  
  virtual matsubara_green_function_t initial_G0(const alps::Parameters& parms, shape_t shape=blockdiagonal);
  matsubara_green_function_t transform_fourpoint(const matsubara_green_function_t& /*G_omega_fourpoint*/,
                                                 const matsubara_green_function_t& G_omega,
                                                 const matsubara_green_function_t& /*G0_omega*/,
                                                 const double /*mu*/, const double /*h*/, const double /*beta*/,
                                                 const unsigned int /*n_matsubara*/) const
  {
    boost::throw_exception(std::logic_error("transform fourpoint function for antiferromagnetic phases not yet implemented"));
    return G_omega;
  }
  
  template<typename T>
  static green_function<T> green_function_inverse(const green_function<T>& G);
  template<typename T>
  static green_function<T> green_function_multiply(const green_function<T>& G1, const green_function<T>& G2);
  template<typename T>
  static green_function<T> green_function_add(const green_function<T>& G1, const green_function<T>& G2);
  template<typename T>
  static green_function<T> green_function_subtract(const green_function<T>& G1, const green_function<T>& G2);
  template<typename T>
  static green_function<T> green_function_divide(const green_function<T>& G1, const T d);
  
  alps::RealObsevaluator kinetic_energy(const std::vector<std::vector<std::vector<alps::RealObsevaluator> > >& selfenergies,
                                        const double beta, const double mu,
                                        const bool error_propagation, const bool paramagnet);
  std::map<double, ClusterHelper::Matrix<std::complex<double> >::Type>
  local_green_function_from_selfenergies(const std::map<double,
                                         std::vector<ClusterHelper::Matrix<std::complex<double> >::Type> >& selfenergies,
                                         const double mu, const double delta);
  
private:
  template <class T>
  green_function<T> spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const;
  template<typename T, class D>
  green_function<T> transform_coordinates(const green_function<T>& G) const;
  
  struct Forward;
  struct Backward;
  
};




//implementation (see also clusterhelper.C, clustertransformer.C)




inline double ClusterHelper::round(const double d) const
{
  int rounded = int(d > 0.0 ? d + 0.5 : d - 0.5);
  if (almost_equal(rounded, d))
    return rounded;
  else
    return d;
}



inline ClusterHelper::vector_type ClusterHelper::round(vector_type v) const
{
  for (unsigned int d=0; d<v.size(); ++d)
    v[d] = round(v[d]);
  return v;
}



inline ClusterHelper::vector_type ClusterHelper::matrix_vector_mult(const Matrix<double>::Type& M, const vector_type& v) const
{
  vector_type vn(v.size());
  for (unsigned int i=0; i<v.size(); ++i) {
    vn[i] = 0.;
    for (unsigned int j=0; j<v.size(); ++j) {
      vn[i] += M(j,i)*v[j];
    }
  }
  return vn;
}


template<typename T>
inline typename ClusterHelper::Matrix<T>::Type ClusterHelper::Matrix<T>::inverse(typename Matrix<T>::Type M)
{
  typename Matrix<T>::Type M_inv(M.size1(), M.size2());
  for (unsigned int i=0; i<M.size1(); ++i)
    for (unsigned int j=0; j<M.size2(); ++j)
      M_inv(i,j) = (T)(i==j);
  boost::numeric::ublas::vector<fortran_int_t> ipivot(M.size1());
  memset(&(ipivot(0)), 0, sizeof(int)*M.size1());
  boost::numeric::bindings::lapack::gesv(M,ipivot, M_inv);
  return M_inv;
}




template<typename T>
inline T ClusterHelper::interpolate(vector_type k, const std::vector<T>& values, const bool antiferromagnet) const
{
  k = round(matrix_vector_mult(reciprocal_basis_inv_, k));
  std::vector<vector_type> neighbor;
  std::vector<unsigned int> neighbor_index;
  switch (dimension()) {
    case 1:
      neighbor.resize(2);
      neighbor_index.resize(2);
      for (int i=0; i<2; ++i)
        neighbor[i].resize(1);
      neighbor[0][0] = floor(k[0]);
      neighbor[1][0] = ceil(k[0]);
      break;
    case 2:
      neighbor.resize(4);
      neighbor_index.resize(4);
      for (int i=0; i<4; ++i)
        neighbor[i].resize(2);
      for (int x=0; x<=1; ++x) {
        double xt = (x==0) ? floor(k[0]) : ceil(k[0]);
        for (int y=0; y<=1; ++y) {
          double yt = (y==0) ? floor(k[1]) : ceil(k[1]);
          neighbor[y*2+x][0] = xt;
          neighbor[y*2+x][1] = yt;
        }
      }
      break;
    case 3:
      neighbor.resize(8);
      neighbor_index.resize(8);
      for (int i=0; i<8; ++i)
        neighbor[i].resize(3);
      for (int x=0; x<=1; ++x) {
        double xt = (x==0) ? floor(k[0]) : ceil(k[0]);
        for (int y=0; y<=1; ++y) {
          double yt = (y==0) ? floor(k[1]) : ceil(k[1]);
          for (int z=0; z<=1; ++z) {
            double zt = (z==0) ? floor(k[2]) : ceil(k[2]);
            neighbor[z*2*2+y*2+x][0] = xt;
            neighbor[z*2*2+y*2+x][1] = yt;
            neighbor[z*2*2+y*2+x][2] = zt;
          }
        }
      }
      break;
  }
  for (std::size_t i=0; i<neighbor.size(); ++i) {
    vector_type nv = matrix_vector_mult(reciprocal_basis_, neighbor[i]);
    nv = map_back_into_brillouin_zone(nv);
    if (antiferromagnet) {
      if (!inside_reduced_brillouin_zone(nv)) {
        for (std::size_t d=0; d<dimension(); ++d)
          nv[d] -= M_PI;
        nv = map_back_into_brillouin_zone(nv);
      }
      if (!inside_reduced_brillouin_zone(nv)) {
        boost::throw_exception(std::logic_error("something doesn't work in ClusterHelper::interpolate)"));
      }
    }
    neighbor_index[i] = identify_vertex(nv, cluster_momenta_);
  }
  std::vector<T> av(neighbor_index.size());
  for (std::size_t i=0; i<neighbor_index.size(); ++i)
    av[i] = values[neighbor_index[i]];
  for (std::size_t d=0; d<dimension(); ++d) {
    double c = k[d] - floor(k[d]) ;
    std::vector<T> t(av.size()/2);
    for (std::size_t i=0; i<av.size(); i+=2) {
      t[i/2] = (1-c)*av[i] + c*av[i+1];
    }
    av = t;
  }
  return av[0];
}



template<typename T>
green_function<T> DCATransformer::symmetrize_helper(const green_function<T>& G,
                                                    const std::pair<std::set<pair_class_type>::const_iterator,
                                                    std::set<pair_class_type>::const_iterator> pair_classes) const
{
  unsigned int nflavor=G.nflavor();
  unsigned int nsite=G.nsite();
  unsigned int nfreq=G.nfreq();
  green_function<T> G_new(nfreq, nsite, nflavor, G.shape());
  std::set<pair_class_type>::const_iterator itc, itc_end;
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    for (spin_t flavor=0; flavor<nflavor; ++flavor) {
      for (frequency_t freq=0; freq<nfreq; ++freq) {
        T av=0;
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          av += G(freq, it->first, it->second, flavor);
        }
        av /= itc->size();
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++)
          G_new(freq, it->first, it->second, flavor) = av;
      }
    }
  }
  return G_new;
}



struct SimpleDCATransformer::Forward {
  template<typename T>
  static void add(green_function<T>& G_new, const green_function<T>& G,
                  const itime_index_t t, const site_t i, const site_t j, const site_t k, const spin_t flavor,
                  double factor) {
    G_new(t, k, k, flavor) += factor * G(t, i, j, flavor);
  }
  template<typename T>
  static void add(nambu_green_function<T>& G_new, const nambu_green_function<T>& G,
                  const itime_index_t t, const site_t i, const site_t j, const site_t k,
                  double factor) {
    G_new(t, k, k, 0,0) += factor * G(t, i, j, 0,0);
    G_new(t, k, k, 0,1) += factor * G(t, i, j, 0,1);
    G_new(t, k, k, 1,0) += factor * G(t, i, j, 1,0);
    G_new(t, k, k, 1,1) += factor * G(t, i, j, 1,1);
  }
};



struct SimpleDCATransformer::Backward {
  template<typename T>
  static void add(green_function<T>& G_new, const green_function<T>& G,
                  const itime_index_t t, const site_t i, const site_t j, const site_t k, const spin_t flavor,
                  double factor) {
    G_new(t, i, j, flavor) += factor * G(t, k, k, flavor);
  }
  template<typename T>
  static void add(nambu_green_function<T>& G_new, const nambu_green_function<T>& G,
                  const itime_index_t t, const site_t i, const site_t j, const site_t k,
                  double factor) {
    G_new(t, i, j, 0, 0) += factor * G(t, k, k, 0, 0);
    G_new(t, i, j, 0, 1) += factor * G(t, k, k, 0, 1);
    G_new(t, i, j, 1, 0) += factor * G(t, k, k, 1, 0);
    G_new(t, i, j, 1, 1) += factor * G(t, k, k, 1, 1);
  }
};



template<typename T, class D>
green_function<T> SimpleDCATransformer::transform_coordinates(const green_function<T>& G) const
{
  unsigned int nflavor=G.nflavor();
  unsigned int nsite=G.nsite();
  unsigned int nt=G.ntime();
  assert(nsite==n_cluster_sites());
  shape_t new_shape=G.shape()==diagonal?nondiagonal:diagonal;
  green_function<T> G_new(nt, nsite, nflavor, new_shape);
  G_new.clear();
  for (site_t k=0; k<nsite; ++k) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
        double phase = 0.;
        for (unsigned int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k,d)*(cluster_coordinate(i,d)-cluster_coordinate(j,d));
        for (itime_index_t t=0; t<nt; ++t) {
          for (spin_t flavor=0; flavor<nflavor; ++flavor) {
            D::add(G_new, G, t, i, j, k, flavor, cos(phase)/double(nsite));
          }
        }
      }
    }
  }
  return G_new;
}


template<class T>
green_function<T> SimpleDCATransformer::spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const
{
  green_function<T> G_new(G);
  if (paramagnet) {
    if(G_new.nflavor()>1 && G_new.nflavor()%2==0) {
      for(spin_t flavor=0; flavor<G_new.nflavor(); flavor+=2) {
        for (site_t k1=0; k1<n_cluster_sites(); ++k1) {
          for (site_t k2=0; k2<n_cluster_sites(); ++k2) {
            for(frequency_t freq=0; freq<G_new.nfreq(); ++freq) {
              G_new(freq, k1, k2, flavor)= 0.5*(G_new(freq, k1, k2, flavor)+G_new(freq, k1, k2, flavor+1));
              G_new(freq, k1, k2, flavor+1)=G_new(freq, k1, k2, flavor);
            }
          }
        }
      }
    }
  }
  return G_new;
}



inline matsubara_green_function_t SimpleDCATransformer::fourpoint_green_function_inverse(const matsubara_green_function_t& G,
                                                                                         const unsigned int nf) const
{
  matsubara_green_function_t G_inv(G);
  for (spin_t f=0; f<G.nflavor(); ++f) {
    Matrix<std::complex<double> >::Type M(nf*G.nsite(), nf*G.nsite());
    for (std::size_t k1=0; k1<G.nsite(); k1++) {
      for (std::size_t k2=0; k2<G.nsite(); k2++) {
        for (std::size_t t1=0; t1<nf; t1++) {
          for (std::size_t t2=0; t2<nf; t2++) {
            M(k1*nf+t1, k2*nf+t2) = G(nf*t1+t2, k1, k2, f);
          }
        }
      }
    }
    /*Matrix<std::complex<double> >::Type M_inv(M.size1(), M.size2());
     for (unsigned int i=0; i<M_inv.size1(); ++i)
     for (unsigned int j=0; j<M_inv.size2(); ++j)
     M_inv(i,j) = (i==j);
     const int n = M.size1();
     const int nrhs = M.size1();
     const int lda = M.size1();
     int ipiv[n];
     const int ldb = M.size1();
     int info;
     zgesv_(&n, &nrhs, (double* )(&M(0,0)), &lda, ipiv, (double* )(&M_inv(0,0)), &ldb, &info);
     std::cout << "(INFO: " << info << ")\n";*/
    Matrix<std::complex<double> >::Type M_inv = Matrix<std::complex<double> >::inverse(M);
    for (std::size_t k1=0; k1<G.nsite(); k1++) {
      for (std::size_t k2=0; k2<G.nsite(); k2++) {
        for (std::size_t t1=0; t1<nf; t1++) {
          for (std::size_t t2=0; t2<nf; t2++) {
            G_inv(nf*t1+t2, k1, k2, f) = M_inv(k1*nf+t1, k2*nf+t2);
          }
        }
      }
    }
  }
  return G_inv;
}



template <typename T>
green_function<T> SimpleDCATransformer::fourpoint_green_function_add(const green_function<T>& G1,
                                                                     const green_function<T>& G2) const
{
  green_function<T> G(G1);
  for (time_t t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (int k1=0; k1<G.nsite(); k1++) {
        for (int k2=0; k2<G.nsite(); k2++) {
          G(t, k1, k2, f) += G2(t, k1, k2, f);
        }
      }
    }
  }
  return G;
}



template <typename T>
green_function<T> SimpleDCATransformer::fourpoint_green_function_subtract(const green_function<T>& G1,
                                                                          const green_function<T>& G2) const
{
  green_function<T> G(G1);
  for (time_t t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (std::size_t k1=0; k1<G.nsite(); k1++) {
        for (std::size_t k2=0; k2<G.nsite(); k2++) {
          G(t, k1, k2, f) -= G2(t, k1, k2, f);
        }
      }
    }
  }
  return G;
}



template <typename T>
green_function<T> SimpleDCATransformer::fourpoint_green_function_divide(const green_function<T>& G1,
                                                                        const T d) const
{
  green_function<T> G(G1);
  for (time_t t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (int k1=0; k1<G.nsite(); k1++) {
        for (int k2=0; k2<G.nsite(); k2++) {
          G(t, k1, k2, f) /= d;
        }
      }
    }
  }
  return G;
}


//////////////AntiFerro DCA////////////////////


template<typename T, class D>
green_function<T> AntiFerroDCATransformer::transform_coordinates(const green_function<T>& G) const
{
  unsigned int nflavor=G.nflavor();
  unsigned int nsite=G.nsite();
  unsigned int nt=G.ntime();
  assert(nsite==n_cluster_sites());
  green_function<T> G_new(nt, nsite, nflavor, blockdiagonal);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
        for (itime_index_t t=0; t<nt; ++t) {
          G_new(t, i, j, flavor) = 0.;
        }
      }
    }
  }
  for (site_t k=0; k<nsite; k+=2) {
    for (site_t q1=0; q1<=1; ++q1) {
      for (site_t q2=0; q2<=1; ++q2) {
        for (site_t i=0; i<nsite; ++i) {
          for (site_t j=0; j<nsite; ++j) {
            double phase = 0.;
            for (unsigned int d=0; d<dimension(); ++d)
              phase += cluster_momentum(k+q1,d)*cluster_coordinate(i,d)-cluster_momentum(k+q2,d)*cluster_coordinate(j,d);
            for (itime_index_t t=0; t<nt; ++t) {
              for (spin_t flavor=0; flavor<nflavor; ++flavor) {
                D::add(G_new, G, t, i, j, k+q1, k+q2, flavor, cos(phase)/double(nsite));
              }
            }
          }
        }
      }
    }
  }
  return G_new;
}


struct AntiFerroDCATransformer::Forward {
  template<typename T>
  static void add(green_function<T>& G_new, const green_function<T>& G,
                  const itime_index_t t, const site_t i, const site_t j, const site_t k1, const site_t k2, const spin_t flavor,
                  double factor)
  {
    G_new(t, k1, k2, flavor) += factor * G(t, i, j, flavor);
  }
};



struct AntiFerroDCATransformer::Backward {
  template<typename T>
  static void add(green_function<T>& G_new, const green_function<T>& G,
                  const itime_index_t t, const site_t i, const site_t j, const site_t k1, const site_t k2, const spin_t flavor,
                  double factor)
  {
    G_new(t, i, j, flavor) += factor * G(t, k1, k2, flavor);
  }
};



template<class T>
green_function<T> AntiFerroDCATransformer::spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const
{
  green_function<T> G_new(G);
  if (paramagnet)
    boost::throw_exception(std::invalid_argument("AFM and PM? redundant parameters! set paramagnet to false.\n"));
  return G_new;
}



template<typename T>
green_function<T> AntiFerroDCATransformer::green_function_inverse(const green_function<T>& G)
{
  green_function<T> G_inv(G);
  for (unsigned int t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (std::size_t k=0; k<G.nsite(); k+=2) {
        T detG = G(t, k, k, f)*G(t, k+1, k+1, f) - G(t, k+1, k, f)*G(t, k, k+1, f);
        G_inv(t, k, k, f) = G(t, k+1, k+1, f)/detG;
        G_inv(t, k+1, k+1, f) = G(t, k, k, f)/detG;
        G_inv(t, k+1, k, f) = -G(t, k+1, k, f)/detG;
        G_inv(t, k, k+1, f) = -G(t, k, k+1, f)/detG;
      }
    }
  }
  return G_inv;
}



template <typename T>
green_function<T> AntiFerroDCATransformer::green_function_multiply(const green_function<T>& G1,
                                                                   const green_function<T>& G2)
{
  green_function<T> G(G1);
  typename Matrix<T>::Type M1(2,2);
  typename Matrix<T>::Type M2(2,2);
  typename Matrix<T>::Type M(2,2);
  for (unsigned int t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (unsigned int k=0; k<G.nsite(); k+=2) {
        M1(0,0) = G1(t, k  , k  , f);
        M1(1,0) = G1(t, k+1, k  , f);
        M1(0,1) = G1(t, k  , k+1, f);
        M1(1,1) = G1(t, k+1, k+1, f);
        M2(0,0) = G2(t, k  , k  , f);
        M2(1,0) = G2(t, k+1, k  , f);
        M2(0,1) = G2(t, k  , k+1, f);
        M2(1,1) = G2(t, k+1, k+1, f);
        M = Matrix<T>::multiply(M1, M2);
        G(t, k, k, f) = M(0,0);
        G(t, k+1, k+1, f) = M(1,1);
        G(t, k+1, k, f) = M(1,0);
        G(t, k, k+1, f) = M(0,1);
      }
    }
  }
  return G;
}



template <typename T>
green_function<T> AntiFerroDCATransformer::green_function_add(const green_function<T>& G1,
                                                              const green_function<T>& G2)
{
  green_function<T> G(G1);
  for (unsigned int t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (unsigned int k=0; k<G.nsite(); k+=2) {
        G(t, k, k, f) += G2(t, k, k, f);
        G(t, k+1, k+1, f) += G2(t, k+1, k+1, f);
        G(t, k+1, k, f) += G2(t, k+1, k, f);
        G(t, k, k+1, f) += G2(t, k, k+1, f);
      }
    }
  }
  return G;
}



template <typename T>
green_function<T> AntiFerroDCATransformer::green_function_subtract(const green_function<T>& G1,
                                                                   const green_function<T>& G2)
{
  green_function<T> G(G1);
  for (unsigned int t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (std::size_t k=0; k<G.nsite(); k+=2) {
        G(t, k, k, f) -= G2(t, k, k, f);
        G(t, k+1, k+1, f) -= G2(t, k+1, k+1, f);
        G(t, k+1, k, f) -= G2(t, k+1, k, f);
        G(t, k, k+1, f) -= G2(t, k, k+1, f);
      }
    }
  }
  return G;
}



template <typename T>
green_function<T> AntiFerroDCATransformer::green_function_divide(const green_function<T>& G1,
                                                                 const T d)
{
  green_function<T> G(G1);
  for (unsigned int t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (unsigned int k=0; k<G.nsite(); k+=2) {
        G(t, k, k, f) /= d;
        G(t, k+1, k+1, f) /= d;
        G(t, k+1, k, f) /= d;
        G(t, k, k+1, f) /= d;
      }
    }
  }
  return G;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////FUNCTIONS FOR NAMBU SPACE//////////////////////////////////////////////////////////////////////////////////


template<typename T, class D>
nambu_green_function<T> SimpleDCATransformer::transform_coordinates(const nambu_green_function<T>& G) const
{
  unsigned int nsite=G.nsite();
  unsigned int nt=G.ntime();
  assert(nsite==n_cluster_sites());
  nambu_green_function<T> G_new(nt, nsite, 2, G.shape()==diagonal?nondiagonal:diagonal);
  G_new.clear();
  //compute G = \sum_{ij} \sum_K exp(i K (R_i-R_j)) G(R_i, R_j)
  for (site_t k=0; k<nsite; ++k) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
        double phase = 0.;
        for (unsigned int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k,d)*(cluster_coordinate(i,d)-cluster_coordinate(j,d));
        for (itime_index_t t=0; t<nt; ++t) {
          D::add(G_new, G, t, i, j, k, cos(phase)/double(nsite));
        }
      }
    }
  }
  return G_new;
}
/*  unsigned int nflavor=G.nflavor();
 unsigned int nsite=G.nsite();
 unsigned int nt=G.ntime();
 assert(nsite==n_cluster_sites());
 green_function<T> G_new(nt, nsite, nflavor);
 G_new.clear();
 for (site_t k=0; k<nsite; ++k) {
 for (site_t i=0; i<nsite; ++i) {
 for (site_t j=0; j<nsite; ++j) {
 double phase = 0.;
 for (unsigned int d=0; d<dimension(); ++d)
 phase += cluster_momentum(k,d)*(cluster_coordinate(i,d)-cluster_coordinate(j,d));
 for (itime_index_t t=0; t<nt; ++t) {
 for (spin_t flavor=0; flavor<nflavor; ++flavor) {
 D::add(G_new, G, t, i, j, k, flavor, cos(phase)/double(nsite));
 }
 }
 }
 }
 }
 return G_new;
 */
template<typename T>
nambu_green_function<T> DCATransformer::symmetrize_helper(const nambu_green_function<T>& G,
                                                          const std::pair<std::set<pair_class_type>::const_iterator,
                                                          std::set<pair_class_type>::const_iterator> pair_classes) const
{
  //  std::cout<<"ignoring symmetrization."<<std::endl;
  //  return G;
  
  unsigned int nsite=G.nsite();
  unsigned int nfreq=G.nfreq();
  nambu_green_function<T> G_new(nfreq, nsite, 2, diagonal);
  G_new=G;
  std::set<pair_class_type>::const_iterator itc, itc_end;
  //std::cout<<"doing spatial symmetries: "<<std::endl;
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    for (spin_t f=0; f<2; ++f) {
      //std::cout<<"new symmetry: "<<std::endl;
      for (frequency_t freq=0; freq<nfreq; ++freq) {
        T av=0;
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          //if(freq==0){ std::cout<<"symmetrizing spatial GF: "<<it->first<<" "<<it->second<<" entry: "<<G(freq, it->first, it->second, f,f)<<std::endl;}
          av += G(freq, it->first, it->second, f,f);
        }
        av /= itc->size();
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++)
          G_new(freq, it->first, it->second, f,f) = av;
      }
    }
  }
  
  //this is the k-space symmetrization respecting d-wave symmetry
  std::vector<double> dwave_prefactor(nsite);
  for(unsigned int i=0;i<nsite;++i){
    double dwave=cos(cluster_momentum(i,0))-cos(cluster_momentum(i,1));
    if(dwave < -1.e-5) dwave_prefactor[i]=-1;
    else if(dwave >  1.e-5) dwave_prefactor[i]=1;
    else dwave_prefactor[i]=0.;
  }
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    //std::cout<<"new symmetry: "<<std::endl;
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      T av01=0;
      T av10=0;
      for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
        if(it->first != it-> second) continue;
        av01+= dwave_prefactor[it->first]*G(freq, it->first, it->second, 0,1);
        av10+= dwave_prefactor[it->first]*G(freq, it->first, it->second, 1,0);
        //if(freq==0){ std::cout<<"symmetrizing spatial GF: "<<it->first<<" "<<it->second<<" dwave: "<<dwave_prefactor[it->first]<<" entry: "<<dwave_prefactor[it->first]*G(freq, it->first, it->second, 0,1)<<" "<<dwave_prefactor[it->first]*G(freq, it->first, it->second, 1,0)<<std::endl;}
      }
      av01/= itc->size();
      av10/= itc->size();
      //if(freq==0){ std::cout<<"setting average to: "<<av01<<" "<<av10<<std::endl; }
      for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
        if(it->first != it-> second) continue;
        G_new(freq, it->first, it->second, 0,1) = dwave_prefactor[it->first]*av01;
        G_new(freq, it->first, it->second, 1,0) = dwave_prefactor[it->first]*av10;
      }
    }
  }
  return G_new;
}


//The CDMFT transormer does the Fourier transforms, dispersions, etc that are needed for the CDMFT.
class CDMFTTransformer : public ClusterTransformer
{
public:
  CDMFTTransformer(const alps::Parameters& params);
  virtual ~CDMFTTransformer(){} //pro forma virtual destructor
  static void generate_transformer(const alps::Parameters& parms, boost::shared_ptr<CDMFTTransformer> &trans_ptr);
  
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
  
  virtual matsubara_green_function_t transform_fourpoint(const matsubara_green_function_t& G_omega_fourpoint,
                                                         const matsubara_green_function_t& G_omega,
                                                         const matsubara_green_function_t& G0_omega,
                                                         const double mu, const double h, const double beta,
                                                         const unsigned int n_matsubara) const =0;
  //Hilbert transform operator
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                                const matsubara_green_function_t &G0_omega,
                                                double &mu, const double h, const double beta)=0;
  //Hilbert transform operator for Nambu GF
  virtual nambu_matsubara_green_function_t operator()(const nambu_matsubara_green_function_t &G_omega,
                                                      const nambu_matsubara_green_function_t &G0_omega,
                                                      nambu_matsubara_green_function_t &selfenergy,
                                                      double &mu, const double h, const double eta1, const double eta2, const double beta, const bool use_selfenergy)=0;
  
  std::complex<double>* exp_iksite(const int sign) const;
  
  int* momentum_translation() const;
  virtual std::vector<double> epsav() const=0;
  virtual std::vector<double> epssqav() const=0;
  double dispersion_moments(const unsigned int i) const { return dispersion_moments_[i]; }
  double bandwidth() const { return bandwidth_; };
  
  double epsilon(const vector_type k) const
  {
    throw std::logic_error("CDMFT dispersion not implemented.");
    return 0.;
  }
  
protected:
  complex_dense_matrix U_;
  
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
  double bandwidth_;
  const int L_;
  double bandwidth_scale_factor_;
  double half_filling_;
  std::vector<double> dispersion_moments_;
};

template<typename T>
green_function<T> CDMFTTransformer::symmetrize_helper(const green_function<T>& G,
                                                      const std::pair<std::set<pair_class_type>::const_iterator,
                                                      std::set<pair_class_type>::const_iterator> pair_classes) const
{
  unsigned int nflavor=G.nflavor();
  unsigned int nsite=G.nsite();
  unsigned int nfreq=G.nfreq();
  green_function<T> G_new(nfreq, nsite, nflavor, nondiagonal);
  std::set<pair_class_type>::const_iterator itc, itc_end;
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    for (spin_t flavor=0; flavor<nflavor; ++flavor) {
      for (frequency_t freq=0; freq<nfreq; ++freq) {
        T av=0;
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++)
          av += G(freq, it->first, it->second, flavor);
        av /= itc->size();
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++)
          G_new(freq, it->first, it->second, flavor) = av;
      }
    }
  }
  return G_new;
}



template<typename T>
nambu_green_function<T> CDMFTTransformer::symmetrize_helper(const nambu_green_function<T>& G,
                                                            const std::pair<std::set<pair_class_type>::const_iterator,
                                                            std::set<pair_class_type>::const_iterator> pair_classes) const
{
  unsigned int nsite=G.nsite();
  unsigned int nfreq=G.nfreq();
  T av;
  nambu_green_function<T> G_new(nfreq, nsite, 2, nondiagonal);
  G_new=G;
  std::set<pair_class_type>::const_iterator itc, itc_end;
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    for (spin_t f=0; f<2; ++f) {
      for (frequency_t freq=0; freq<nfreq; ++freq) {
        av=0;
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          av += G(freq, it->first, it->second, f,f);
        }
        av /= itc->size();
        for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
          G_new(freq, it->first, it->second, f,f) = av;
        }
      }
    }
    std::cout<<std::endl;
  }
  
  //transform the dwave factor of K into realspace
  dense_matrix dwave(nsite,nsite);dwave.clear();
  dwave(0,1)=-1;dwave(0,3)=1;dwave(1,0)=-1;dwave(1,2)=1;dwave(2,1)=1;dwave(2,3)=-1;dwave(3,0)=1;dwave(3,2)=-1;
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    //std::cout<<"new symmetry: "<<std::endl;
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      T av01=0;
      T av10=0;
      for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
        //if(it->first != it-> second) continue;
        av01+= dwave(it->first,it->second)*G(freq, it->first, it->second, 0,1);
        av10+= dwave(it->first,it->second)*G(freq, it->first, it->second, 1,0);
        //if(freq==0){ std::cout<<"symmetrizing spatial GF: "<<it->first<<" "<<it->second<<" dwave: "<<dwave_prefactor[it->first]<<" entry: "<<dwave_prefactor[it->first]*G(freq, it->first, it->second, 0,1)<<" "<<dwave_prefactor[it->first]*G(freq, it->first, it->second, 1,0)<<std::endl;}
      }
      av01/= itc->size();
      av10/= itc->size();
      //if(freq==0){ std::cout<<"setting average to: "<<av01<<" "<<av10<<std::endl; }
      for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++){
        //if(it->first != it-> second) continue;
        G_new(freq, it->first, it->second, 0,1) = dwave(it->first,it->second)*av01;
        G_new(freq, it->first, it->second, 1,0) = dwave(it->first,it->second)*av10;
      }
    }
  }
  return G_new;
}


//does the CDMFT cluster transform for the paramagnetic phase (no AFM order)
class FoursiteCDMFTTransformer : public CDMFTTransformer
{
public:
  FoursiteCDMFTTransformer(const alps::Parameters& params) :
  CDMFTTransformer(params),
  compute_filling_(params.defined("FILLING")),
  n_target_(compute_filling_?(double)(params["FILLING"]):-1),
  epsav_(n_cluster_sites(),0.),
  epssqav_(n_cluster_sites(),0.),
  L_(params["L"]),
  t_(params["t"]),
  tprime_(params.value_or_default("tprime",0)),
  tstar_(params["tstar"])
  {
    params_ptr_=&params;
    set_U_(n_cluster_sites(),U_);
    /*  U_.resize(4,4);
     U_(0,0)=1.; U_(0,1)=1.; U_(0,2)=1.; U_(0,3)=1.;
     U_(1,0)=1.; U_(1,1)=-1.; U_(1,2)=1.; U_(1,3)=-1.;
     U_(2,0)=1.; U_(2,1)=1.; U_(2,2)=-1.; U_(2,3)=-1.;
     U_(3,0)=1.; U_(3,1)=-1.; U_(3,2)=-1.; U_(3,3)=1.;
     U_/=2.;
     
     /*epsav_[0]=-8*(M_PI*t+2*tprime)/(M_PI*M_PI);
     epsav_[1]= 8*(M_PI*t-2*tprime)/(M_PI*M_PI);
     epsav_[2]=16*tprime/(M_PI*M_PI);
     epsav_[3]=16*tprime/(M_PI*M_PI);
     
     epssqav_[0]=4*((8+M_PI*M_PI)*t*t+8*M_PI*t*tprime+M_PI*M_PI*tprime*tprime)/(M_PI*M_PI);
     epssqav_[1]=4*((8+M_PI*M_PI)*t*t-8*M_PI*t*tprime+M_PI*M_PI*tprime*tprime)/(M_PI*M_PI);
     epssqav_[2]=4*((1-8/(M_PI*M_PI))*t*t+tprime*tprime);
     epssqav_[3]=4*((1-8/(M_PI*M_PI))*t*t+tprime*tprime);*/
  }
  //Hilbert transform operator
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
  //Fourier transform functions
  matsubara_green_function_t transform_into_real_space(const matsubara_green_function_t& G) const;
  matsubara_green_function_t transform_into_k_space(const matsubara_green_function_t& G) const ;
  itime_green_function_t transform_into_real_space(const itime_green_function_t& G) const ;
  itime_green_function_t transform_into_k_space(const itime_green_function_t& G) const ;
  
  //spin symmetrization functions in Nambu space
  nambu_itime_green_function_t symmetrize(const nambu_itime_green_function_t& G, const bool paramagnet) const
  {
    return spin_symmetrizer(G, paramagnet);
  }
  nambu_matsubara_green_function_t symmetrize(const nambu_matsubara_green_function_t& G, const bool paramagnet) const
  {
    return spin_symmetrizer(G, paramagnet);
  }
  nambu_matsubara_green_function_t transform_into_real_space(const nambu_matsubara_green_function_t& G) const ;
  nambu_matsubara_green_function_t transform_into_k_space(const nambu_matsubara_green_function_t& G) const;
  nambu_itime_green_function_t transform_into_real_space(const nambu_itime_green_function_t& G) const ;
  nambu_itime_green_function_t transform_into_k_space(const nambu_itime_green_function_t& G) const ;
  
  //initial G0
  virtual matsubara_green_function_t initial_G0      (const alps::Parameters& parms, shape_t shape=nondiagonal);
  virtual nambu_matsubara_green_function_t initial_G0_nambu(const alps::Parameters& parms);
  matsubara_green_function_t transform_fourpoint(const matsubara_green_function_t& G_omega_foupoint,
                                                 const matsubara_green_function_t &G_omega,
                                                 const matsubara_green_function_t &G0_omega,
                                                 const double mu, const double h, const double beta,
                                                 const unsigned int n_matsubara) const;
  std::map<double, std::complex<double> > local_green_function_from_selfenergies(const std::map<double, std::vector<std::complex<double> > >&, double, double);
  alps::RealObsevaluator kinetic_energy(const std::vector<std::vector<std::vector<alps::SimpleObservableEvaluator<double> > > >&, double, double, bool, bool);
  std::vector<double> epsav() const;
  std::vector<double> epssqav() const;
private:
  template<typename T, class D>
  green_function<T> transform_coordinates(const       green_function<T>& G) const;
  template<typename T, class D>
  nambu_green_function<T> transform_coordinates(const nambu_green_function<T>& G) const;
  template <class T>
  green_function<T> spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const;
  
  //in nambu case: <double> very different from <complex<double> > -> no template
  nambu_green_function<std::complex<double> > spin_symmetrizer(const nambu_green_function<std::complex<double> >& G, const bool paramagnet) const{ throw std::runtime_error("not implemented.");}
  nambu_green_function<double> spin_symmetrizer(const nambu_green_function<double>& G, const bool paramagnet) const{ throw std::runtime_error("not implemented.");}
  
  matsubara_green_function_t fourpoint_green_function_inverse(const matsubara_green_function_t& G, const unsigned int nt) const;
  template<typename T>
  green_function<T> fourpoint_green_function_add(const green_function<T>& G1,
                                                 const green_function<T>& G2) const;
  template<typename T>
  green_function<T> fourpoint_green_function_subtract(const green_function<T>& G1,
                                                      const green_function<T>& G2) const;
  template<typename T>
  green_function<T> fourpoint_green_function_divide(const green_function<T>& G1,
                                                    const T d) const;
  void set_Tij_matrix(double kx, double ky, complex_dense_matrix &H0);
  void set_nambu_Tij_dwave_matrix(double kx, double ky, double eta1,double eta2,double mu, complex_dense_matrix &H0);
  void set_U_(unsigned int nsite, complex_dense_matrix &U_);
  struct Forward;
  struct Backward;
  
  bool compute_filling_;
  double n_target_;
  const alps::Parameters *params_ptr_;
  //complex_dense_matrix U_;
  std::vector<double> epsav_;
  std::vector<double> epssqav_;
  const int L_;
  const double t_;
  const double tprime_;
  const double tstar_;
};

template<class T>
green_function<T> FoursiteCDMFTTransformer::spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const
{
  green_function<T> G_new(G);
  if (paramagnet) {
    if(G_new.nflavor()>1 && G_new.nflavor()%2==0) {
      for(spin_t flavor=0; flavor<G_new.nflavor(); flavor+=2) {
        for (site_t k1=0; k1<n_cluster_sites(); ++k1) {
          for (site_t k2=0; k2<n_cluster_sites(); ++k2) {
            for(frequency_t freq=0; freq<G_new.nfreq(); ++freq) {
              G_new(freq, k1, k2, flavor)= 0.5*(G_new(freq, k1, k2, flavor)+G_new(freq, k1, k2, flavor+1));
              G_new(freq, k1, k2, flavor+1)=G_new(freq, k1, k2, flavor);
            }
          }
        }
      }
    }
  }
  return G_new;
}

#endif /*ALPS_CLUSTER_H*/
