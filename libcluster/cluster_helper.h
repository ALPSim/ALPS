/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2007      by Sebastian Fuchs <fuchs@comp-phys.org>
 * Copyright (C) 2007-2015 by Emanuel Gull <egull@umich.edu>
 *
 *
 *****************************************************************************/

#pragma once

#include <set>
#include "boost/tuple/tuple.hpp"

#include "alps/params.hpp"
#include "alps/hdf5/pointer.hpp"

#include "blasheader.h"
#include "general_matrix.h"
#include "dump.h"

#include "cluster_graph.h"


class ClusterHelper : private alps_graph
{
public:
  
  template <typename T> struct Matrix
  {
    typedef blas::general_matrix<T> Type;
    
    static Type inverse(Type M){return M.invert();}
    
    static Type multiply(Type M1, Type M2) {
      return M1*M2;
    }
    static double norm_1(const Type &M){ return M.max();}
  };
  static void define_parameters(alps::params &p);
  
  typedef std::set<int> site_class_type;
  typedef std::pair<int, int> pair_type;
  typedef std::set<pair_type> pair_class_type;
  typedef std::set<site_class_type> site_class_set_type;
  typedef std::set<pair_class_type> pair_class_set_type;
  typedef std::pair<site_class_set_type::const_iterator, site_class_set_type::const_iterator> site_class_set_itpair_type;
  typedef std::pair<pair_class_set_type::const_iterator, pair_class_set_type::const_iterator> pair_class_set_itpair_type;
  
  ClusterHelper(const alps::params& params);
  virtual ~ClusterHelper(){}
  
  double cluster_coordinate(int i, int j) const { return cluster_coordinates_(i,j); }
  double cluster_momentum(int i, int j) const { return cluster_momenta_(i,j); }
  double lattice_momentum( int i,  int j) const { return lattice_momenta_(i,j); }
  double reciprocal_basis_vector(const unsigned int i, const unsigned int j) const { return reciprocal_basis_(i,j); }
  double basis_vector(const unsigned int i, const unsigned int j) const { return basis_(i,j); }
  ///getter function for periodicity (number of hops in each direction until we're back on the same site)
  int periodicity(const unsigned int i) const { return periodicity_[i]; }
  bool bipartite() const { return bipartite_; }
  int n_cluster_sites() const { return alps_graph::unit_cell_volume(); }
  int n_cluster_cells() const { return alps_graph::volume(); }
  int n_lattice_momenta() const { return lattice_momenta_.size1(); }
  int dimension() const { return alps_graph::dimension(); }
  int cluster_momentum_approximation(std::vector<double> k) const;
  
  template<typename T>
  T interpolate(std::vector<double> k, const std::vector<T>& values, const bool antiferromagnet) const;
  double site_distance(const unsigned int i, const unsigned int j) const;
  
  site_class_set_itpair_type site_classes_real_space() const
  { return site_class_set_itpair_type (site_class_set_real_space_.begin(), site_class_set_real_space_.end()); }
  pair_class_set_itpair_type pair_classes_real_space() const
  { return pair_class_set_itpair_type(pair_class_set_real_space_.begin(), pair_class_set_real_space_.end()); }
  site_class_set_itpair_type site_classes_k_space() const
  { return site_class_set_itpair_type(site_class_set_k_space_.begin(), site_class_set_k_space_.end()); }
  pair_class_set_itpair_type pair_classes_k_space() const
  { return pair_class_set_itpair_type(pair_class_set_k_space_.begin(), pair_class_set_k_space_.end()); }
  const Matrix<int>::Type &symmetry_table_k_space() const{return symmetry_table_k_space_;}
  const pair_class_set_type &pair_class_set_k_space() const{return pair_class_set_k_space_;}
  
  virtual void write_hdf5(alps::hdf5::archive &ar, const std::string &p) const;
  bool write() const;
  bool read();
  
  bool operator==(const ClusterHelper &CH) const;
private:
  void setup_tables(const alps::params& params);
  ///create a table of point group symmetries in 1d, 2d, or 3d.
  std::vector<Matrix<double>::Type > setup_point_groups() const;
  void symmetry_analysis(const Matrix<double>::Type& basis,
                         const Matrix<double>::Type& basis_inv);
//  void ladder_symmetry_analysis();
  std::pair<site_class_set_type, pair_class_set_type> setup_classes(const Matrix<int>::Type& symmetry_table) const;
  bool inside_cluster(std::vector<double> v, const Matrix<double>::Type& b) const;
  bool inside_reduced_brillouin_zone(std::vector<double> v) const;
  bool append_column_no_duplicates(Matrix<double>::Type& m, const std::vector<double>& v) const;
  bool append_column_no_duplicates(Matrix<int>::Type& m, std::map<int, int>& smap) const;
  std::vector<double> map_back_into_brillouin_zone(std::vector<double> v) const;
  std::vector<double> map_back_into_cluster(std::vector<double> v,
                                    const Matrix<double>::Type& B_inv, const Matrix<double>::Type& B) const;
  unsigned int identify_vertex(const std::vector<double>& v, const Matrix<double>::Type& coordinate_matrix) const;
  unsigned int identify_sublattice(const Matrix<double>::Type& coordinate_matrix, const unsigned int i) const;
  double round(const double d) const;
  std::vector<double> round(std::vector<double> v) const;
  bool almost_equal(const double d1, const double d2) const { return std::abs(d1-d2)<epsilon_; }
  std::vector<double> matrix_vector_mult(const Matrix<double>::Type& M, const std::vector<double>& v) const;
  void write_site_class_set(alps::hdf5::archive &ar, const site_class_set_type &set, const std::string &name,const std::string &location="/transformer") const;
  void write_pair_class_set(alps::hdf5::archive &ar, const pair_class_set_type &set, const std::string &name,const std::string &location="/transformer") const;
  void read_site_class_set(alps::hdf5::archive &ar, site_class_set_type &set, const std::string &name,const std::string &location="/transformer");
  void read_pair_class_set(alps::hdf5::archive &ar, pair_class_set_type &set, const std::string &name,const std::string &location="/transformer");
  void find_non_bipartite_cluster_momenta();
  void find_bipartite_cluster_momenta();
  void find_integration_weights();
  ///Go through the cluster in each direction and find how many hops it takes to get to the same place
  void find_periodicity();
  ///detect if the cluster is bipartite or not
  void detect_if_bipartite();
  ///compute the cluter coordinates and store them in the 'cluster_coordinates' matrix.
  void extract_and_arrange_cluster_coordinates();
  ///transcribe the basis form alps_graph and invert it
  void compute_inverse_basis();
  ///transcribe reciprocal basis form alps_graph and invert it
  void compute_inverse_reciprocal_basis();


  Matrix<double>::Type cluster_coordinates_;
  Matrix<double>::Type cluster_momenta_;
  Matrix<double>::Type lattice_momenta_;
  Matrix<double>::Type reciprocal_basis_;
  Matrix<double>::Type reciprocal_basis_inv_;
  Matrix<double>::Type rbzone_boundaries_;
  Matrix<double>::Type basis_;
  Matrix<double>::Type basis_inv_;
  site_class_set_type site_class_set_real_space_;
  pair_class_set_type pair_class_set_real_space_;
  site_class_set_type site_class_set_k_space_;
  pair_class_set_type pair_class_set_k_space_;
  Matrix<int>::Type symmetry_table_k_space_;
  Matrix<int>::Type symmetry_table_real_space_;
  ///number of hops in each direction until we're back at the same place
  std::vector<int> periodicity_;

  bool paramagnet_;
  bool bipartite_;
  const double epsilon_;
  const std::string hdf5_dump_name_;
  const std::string cluster_name_;
};


