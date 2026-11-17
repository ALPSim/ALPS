/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2007 by Sebastian Fuchs <fuchs@comp-phys.org>
 * Copyright (C) 2007-2015 by Emanuel Gull <egull@umich.edu>
 *
 *
 *****************************************************************************/

#include <alps/config.hpp> // to know if we have MPI

#ifdef ALPS_HAVE_MPI
#include <mpi.h>
#endif

#include "cluster_graph.h"
#include "cluster_helper.h"

void ClusterHelper::define_parameters(alps::params &p){
  p.define<std::string>("dca.LATTICE", "name of the cluster");
  p.define<std::string>("dca.LATTICE_LIBRARY", "name of the xml cluster definition file");
  p.define<std::string>("dca.LATTICE_DUMP", "cluster.h5", "name of the cluster transformer hdf5 dump file");
  p.define<bool>("dca.AFM_SELFCONSISTENCY", false, "Use the AFM double unit cell framework");
  p.define<bool>("dca.CO_SELFCONSISTENCY", false, "Use the charge order double unit cell framework.");
  p.define<int>("dca.SITES",  "number of cluster sites");
  p.define<int>("dca.L", "linear discretization of cluster k-point integration");
}
ClusterHelper::ClusterHelper(const alps::params& params) :
alps_graph(params["dca.LATTICE"],params["dca.LATTICE_LIBRARY"],params["dca.L"]),
paramagnet_(!params["dca.AFM_SELFCONSISTENCY"]),
epsilon_(1e-9),
hdf5_dump_name_(params["dca.LATTICE_DUMP"].as<std::string>()),
cluster_name_(params["dca.LATTICE"].as<std::string>())
{
  //check that the numer of sites is consistent between the 'SITES' parameter and what we get from the cluster.
  if (int(params["dca.SITES"]) != alps_graph::unit_cell_volume()){
   std::stringstream error_message; error_message<<"the graph: "<<params["dca.LATTICE"].as<std::string>()<<" unit cell volume is not :"<<params["dca.SITES"].as<int>();
   throw std::runtime_error(error_message.str());
  }
  //generate the symmetry tables for this transformer.
  if(!read())
    setup_tables(params);
  
  //Find the periodicity in each direction
  find_periodicity();
}

//Find the periodicity (i.e. the number of hops until we're back on the same site) in each direction
void ClusterHelper::find_periodicity() {
  periodicity_.clear();
  periodicity_.resize(dimension(), 0.);
  for (int k = 0; k < dimension(); k++) {
    int n = 0;
    bool done = false;
    while (!done) {
      n += 1;
      for (int d = 0; d < dimension(); d++) {
        if (n * reciprocal_basis_(k, d) > 2 * M_PI) {
          done = true;
        }
      }
    }
    periodicity_[k] = n - 1;
  }
}

void ClusterHelper::find_non_bipartite_cluster_momenta() {
  for(int orig=0; orig<cell_origins().size(); ++orig){
    std::vector<double> o = cell_origins()[orig];
    for (unsigned int i = 0; i < n_cluster_sites(); ++i) {
      std::vector<double> c = o;
      for (unsigned int d = 0; d < dimension(); ++d)
        c[d] += cluster_coordinate(i, d);
      c = matrix_vector_mult(reciprocal_basis_, c);
      c = map_back_into_brillouin_zone(c);
      append_column_no_duplicates(cluster_momenta_, c);
    }
  }
}

void ClusterHelper::find_bipartite_cluster_momenta() {
  switch (dimension()) {
    case 1:
      rbzone_boundaries_.resize(1, 1);
      rbzone_boundaries_(0, 0) = M_PI / 2;
      break;
    case 2:
      rbzone_boundaries_.resize(2, 2);
      rbzone_boundaries_(0, 0) = M_PI / 2;
      rbzone_boundaries_(0, 1) = M_PI / 2;
      rbzone_boundaries_(1, 0) = -M_PI / 2;
      rbzone_boundaries_(1, 1) = M_PI / 2;
      break;
    case 3:
      //reduced Brillouin zone is the one of an fcc lattice with a=2
      rbzone_boundaries_.resize(4, 3);
      rbzone_boundaries_(0, 0) = M_PI / 2;
      rbzone_boundaries_(0, 1) = M_PI / 2;
      rbzone_boundaries_(0, 2) = M_PI / 2;
      rbzone_boundaries_(1, 0) = -M_PI / 2;
      rbzone_boundaries_(1, 1) = M_PI / 2;
      rbzone_boundaries_(1, 2) = M_PI / 2;
      rbzone_boundaries_(2, 0) = M_PI / 2;
      rbzone_boundaries_(2, 1) = -M_PI / 2;
      rbzone_boundaries_(2, 2) = M_PI / 2;
      rbzone_boundaries_(3, 0) = M_PI / 2;
      rbzone_boundaries_(3, 1) = M_PI / 2;
      rbzone_boundaries_(3, 2) = -M_PI / 2;
      break;
  }
  for(int orig=0; orig<cell_origins().size(); ++orig){
    std::vector<double> o = cell_origins()[orig];
    for (unsigned int i = 0; i < n_cluster_sites(); ++i) {
      std::vector<double> c = o;
      for (unsigned int d = 0; d < dimension(); ++d)
        c[d] += cluster_coordinate(i, d);
      c = matrix_vector_mult(reciprocal_basis_, c);
      c = map_back_into_brillouin_zone(c);
      if (inside_reduced_brillouin_zone(c)) {
        append_column_no_duplicates(cluster_momenta_, c);
        for (std::size_t d = 0; d < dimension(); ++d)
          c[d] += M_PI;
        c = map_back_into_brillouin_zone(c);
        append_column_no_duplicates(cluster_momenta_, c);
      }
    }
  }
}

void ClusterHelper::find_integration_weights() {
  //find lattice integration weights (Simpson weights)
  std::vector<double> weights;
  for(int i=0;i<cell_momenta().size();++i){
    //for (; itm != itm_end; ++itm) {
    std::vector<double> mom = cell_momenta()[i];
    std::vector<double> mom_trans = matrix_vector_mult(reciprocal_basis_inv_, mom);
    std::vector<std::vector<double> > moms;
    switch (dimension()) {
      case 1:
        for (int s1 = -1; s1 <= 1; s1 += 2) {
          std::vector<double> tmom = mom_trans;
          tmom[0] = s1 * mom_trans[0];
          moms.push_back(tmom);
        }
        break;
      case 2:
        for (int s1 = -1; s1 <= 1; s1 += 2) {
          for (int s2 = -1; s2 <= 1; s2 += 2) {
            std::vector<double> tmom = mom_trans;
            tmom[0] = s1 * mom_trans[0];
            tmom[1] = s2 * mom_trans[1];
            moms.push_back(tmom);
          }
        }
        break;
      case 3:
        for (int s1 = -1; s1 <= 1; s1 += 2) {
          for (int s2 = -1; s2 <= 1; s2 += 2) {
            for (int s3 = -1; s3 <= 1; s3 += 2) {
              std::vector<double> tmom = mom_trans;
              tmom[0] = s1 * mom_trans[0];
              tmom[1] = s2 * mom_trans[1];
              tmom[2] = s3 * mom_trans[2];
              moms.push_back(tmom);
            }
          }
        }
        break;
    }
    double L = pow(n_cluster_cells(), 1. / double(dimension()));
    for (std::size_t g = 0; g < moms.size(); ++g) {
      double weight = 1.;
      for (std::size_t d = 0; d < dimension(); ++d) {
        double m = (moms[g][d] + 0.5) * L;
        if (almost_equal(m, 0) || almost_equal(m, L)) {
          //boundary gets a third
          weight *= 1. / 3.;
        } else
          //rest gets 2/3 and 4/3 of weight according to Simpson
          if ((((int) (std::abs(round(m)))) % 2 == 0
               && ((int) (std::abs(round(L)))) % 2 == 0)
              || (((int) (std::abs(round(m)))) % 2 == 1
                  && ((int) (std::abs(round(L)))) % 2 == 1)) {
                weight *= 2. / 3.;
              } else {
                weight *= 4. / 3.;
              }
      }
      //momw has dimension+1 entries: dimension for the actual coordinate, and the last one for the weight
      std::vector<double> momw(dimension() + 1);
      std::vector<double> momt = matrix_vector_mult(reciprocal_basis_, moms[g]);
      for (std::size_t d = 0; d < dimension(); ++d)
        momw[d] = momt[d];
      momw[dimension()] = weight;
      //std::cout<<"appending momentum to lattice_momenta_, weight: "<<weight<<std::endl;
      append_column_no_duplicates(lattice_momenta_, momw);
    }
  }
}

void ClusterHelper::detect_if_bipartite() {
  int sum = 0;
  for (unsigned int i = 0; i < basis_vectors().size(); ++i) {
    int sum_i = 0;
    for (unsigned int d = 0; d < dimension(); ++d) {
      sum_i += int(round(basis_vectors()[i][d]));
    }
    sum += sum_i >= 0 ? sum_i % 2 : (-sum_i) % 2;
  }
  bipartite_ = (sum == 0);
}
void ClusterHelper::compute_inverse_basis(){
  basis_.resize(dimension(), dimension());
  for (unsigned int i=0; i<basis_vectors().size(); ++i) {
    for (unsigned int d=0; d<dimension(); ++d) {
      basis_(i,d) = basis_vectors()[i][d];
    }
  }
  basis_inv_ = Matrix<double>::inverse(basis_);
}
void ClusterHelper::compute_inverse_reciprocal_basis(){
  reciprocal_basis_.resize(dimension(), dimension());
  for (unsigned int i=0; i<reciprocal_basis_vectors().size();++i){
    for (unsigned int d=0; d<dimension(); ++d){
      reciprocal_basis_(i,d) = reciprocal_basis_vectors()[i][d];
    }
  }
  reciprocal_basis_inv_ = Matrix<double>::inverse(reciprocal_basis_);
}
void ClusterHelper::extract_and_arrange_cluster_coordinates(){
  int maxc=n_cluster_sites();
  std::vector<double> c(dimension());
  switch(dimension()) {
    case 1 :
      for (int i=-maxc; i<=maxc; ++i) {
        c[0] = i;
        if (inside_cluster(c, basis_inv_))
          append_column_no_duplicates(cluster_coordinates_, c);
      }
      break;
    case 2:
      for (int i=-maxc; i<=maxc; ++i) {
        for (int j=-maxc; j<=maxc; ++j) {
          c[0] = i; c[1] = j;
          if (inside_cluster(c, basis_inv_))
            append_column_no_duplicates(cluster_coordinates_, c);
        }
      }
      break;
    case 3:
      for (int i=-maxc; i<=maxc; ++i) {
        for (int j=-maxc; j<=maxc; ++j) {
          for (int k=-maxc; k<=maxc; ++k) {
            c[0] = i; c[1] = j; c[2]=k;
            if (inside_cluster(c, basis_inv_)) {
              append_column_no_duplicates(cluster_coordinates_, c);
            }
          }
        }
      }
      break;
  }

  //arrange vertices for bipartite lattices
  if (bipartite()) {
    for (std::size_t i=0; i<cluster_coordinates_.size1(); ++i) {
      if (identify_sublattice(cluster_coordinates_, i)!=i%2) {
        std::size_t j = i;
        while (identify_sublattice(cluster_coordinates_, ++j)!=i%2) {};
        for (std::size_t d=0; d<cluster_coordinates_.size2(); ++d) {
          double t=cluster_coordinates_(i,d);
          cluster_coordinates_(i,d) = cluster_coordinates_(j,d);
          cluster_coordinates_(j,d) = t;
        }
      }
    }
  }

}
void ClusterHelper::setup_tables(const alps::params &params)
{
  detect_if_bipartite();
  
  compute_inverse_basis();
  
  extract_and_arrange_cluster_coordinates();
  
  compute_inverse_reciprocal_basis();
  
  //setup cluster momenta
  if (bipartite()) {
    find_bipartite_cluster_momenta();
  }
  else {
    find_non_bipartite_cluster_momenta();
  }
  
  //find lattice integration weights (Simpson weights)
  find_integration_weights();
  
  if (cluster_coordinates_.size1() != n_cluster_sites())
    boost::throw_exception(std::logic_error("parameter SITES does not match cluster definition"));
  if (cluster_momenta_.size1() != n_cluster_sites())
    boost::throw_exception(std::logic_error("setup of cluster momenta table failed"));
  symmetry_analysis(basis_, basis_inv_);
#ifndef ALPS_HAVE_MPI
  write();
#else
  int mpi_init; MPI_Initialized(&mpi_init);
  if(mpi_init){
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0){
    write();
  }
  }else{write();}
#endif
}
void ClusterHelper::write_site_class_set(alps::hdf5::archive &ar, const site_class_set_type &set, const std::string &name, const std::string &location) const{
  int z=0;
  std::stringstream path; path<<"/"+location+"/"<<name<<"/size";
  ar<<alps::make_pvp(path.str(),set.size());
  for(site_class_set_type::const_iterator it=set.begin(); it != set.end();++it){
    int z1=0;
    std::stringstream path; path<<"/"+location+"/"<<name<<"/"<<z<<"/size";
    ar<<alps::make_pvp(path.str(),it->size());
    for(site_class_type::const_iterator it1=it->begin(); it1 != it->end();++it1){
      std::stringstream path; path<<"/"+location+"/"<<name<<"/"<<z<<"/"<<z1;
      ar<<alps::make_pvp(path.str(), *it1);
      z1++;
    }
    z++;
  }
}
void ClusterHelper::write_pair_class_set(alps::hdf5::archive &ar, const pair_class_set_type &set, const std::string &name, const std::string &location) const{
  int z1=0;
  std::vector<int> indices;
  std::vector<int> values;
  for(pair_class_set_type::const_iterator it=set.begin(); it != set.end();++it){
    indices.push_back(z1);
    for(pair_class_type::const_iterator it1=it->begin(); it1 != it->end();++it1){
      values.push_back(it1->first);
      values.push_back(it1->second);
      z1++;
    }
  }
  indices.push_back(z1);
  std::stringstream path1, path2;
  path1<<"/"+location+"/"<<name<<"/indices";
  ar<<alps::make_pvp(path1.str(), indices);
  path2<<"/"+location+"/"<<name<<"/values";
  ar<<alps::make_pvp(path2.str(), values);
}
void ClusterHelper::read_site_class_set(alps::hdf5::archive &ar, site_class_set_type &set, const std::string &name, const std::string &location){
  int zmax, z1max;
  std::stringstream path; path<<"/"+location+"/"<<name<<"/size";
  ar>>alps::make_pvp(path.str(),zmax);
  for(int z=0;z<zmax;++z){
    std::stringstream path; path<<"/"+location+"/"<<name<<"/"<<z<<"/size";
    ar>>alps::make_pvp(path.str(),z1max);
    std::set<int> s;
    for(int z1=0;z1<z1max;++z1){
      int val;
      std::stringstream path; path<<"/"+location+"/"<<name<<"/"<<z<<"/"<<z1;
      ar>>alps::make_pvp(path.str(), val);
      s.insert(val);
    }
    set.insert(s);
  }
}
void ClusterHelper::read_pair_class_set(alps::hdf5::archive &ar, pair_class_set_type &set, const std::string &name, const std::string &location){
  std::stringstream path; path<<"/"+location+"/"<<name;
  std::string root = path.str();
  std::vector<int> indices;
  ar>>alps::make_pvp(root + "/indices",indices);
  std::vector<int> values;
  ar>>alps::make_pvp(root + "/values",values);
  for(int z = 0; z < indices.size()-1; ++z) {
    std::set<pair_type> s;
    for(int z1 = indices[z]; z1 < indices[z+1]; ++z1) {
      pair_type val;
      val.first  = values[z1*2    ];
      val.second = values[z1*2 + 1];
      s.insert(val);
    }
    set.insert(s);
  }
}

bool ClusterHelper::write() const{
  std::cout<<"entering write for name: "<<hdf5_dump_name_<<" cluster name: "<<cluster_name_<<std::endl;
    alps::hdf5::archive ar(hdf5_dump_name_, "w");
    if(ar.is_data("/"+cluster_name_+"/name")) return false; //already have a transformer there.
    ar<<alps::make_pvp("/"+cluster_name_+"/name", cluster_name_);
    ar<<alps::make_pvp("/"+cluster_name_+"/bipartite", bipartite_);
    ar<<alps::make_pvp("/"+cluster_name_+"/epsilon", epsilon_);
    ar<<alps::make_pvp("/"+cluster_name_+"/cluster_coordinates", cluster_coordinates_);
    ar<<alps::make_pvp("/"+cluster_name_+"/cluster_momenta", cluster_momenta_);
    ar<<alps::make_pvp("/"+cluster_name_+"/lattice_momenta", lattice_momenta_);
    ar<<alps::make_pvp("/"+cluster_name_+"/reciprocal_basis", reciprocal_basis_);
    ar<<alps::make_pvp("/"+cluster_name_+"/reciprocal_basis_inv", reciprocal_basis_inv_);
    ar<<alps::make_pvp("/"+cluster_name_+"/rbzone_boundaries", rbzone_boundaries_);
    ar<<alps::make_pvp("/"+cluster_name_+"/symmetry_table_k_space", symmetry_table_k_space_);
    ar<<alps::make_pvp("/"+cluster_name_+"/symmetry_table_real_space", symmetry_table_real_space_);
    
    //the sets we need to do by hand:
    //site class set real space
    write_site_class_set(ar, site_class_set_real_space_, "site_class_set_real_space", cluster_name_);
    write_site_class_set(ar, site_class_set_k_space_, "site_class_set_k_space", cluster_name_);
    write_pair_class_set(ar, pair_class_set_real_space_, "pair_class_set_real_space", cluster_name_);
    write_pair_class_set(ar, pair_class_set_k_space_, "pair_class_set_k_space", cluster_name_);
    std::cout<<"exiting write for name: "<<hdf5_dump_name_<<" cluster name: "<<cluster_name_<<std::endl;
    return true;
}
bool ClusterHelper::read(){
    //test if we can read the file
    {
      std::ifstream in_file(hdf5_dump_name_.c_str());
      if(!in_file.is_open()) return false;
    }
    alps::hdf5::archive ar(hdf5_dump_name_, "r");
    if(!ar.is_data("/"+cluster_name_+"/name")) return false;
    std::string name;     ar>>alps::make_pvp("/"+cluster_name_+"/name", name);
    if(name!=cluster_name_) throw std::logic_error("reading cluster transformer with name: "+name+" for cluster: "+cluster_name_);
    //write data to the 'tmpname' file
    ar>>alps::make_pvp("/"+cluster_name_+"/bipartite", bipartite_);
    ar>>alps::make_pvp("/"+cluster_name_+"/epsilon", *(const_cast<double*>(&epsilon_)));
    ar>>alps::make_pvp("/"+cluster_name_+"/cluster_coordinates", cluster_coordinates_);
    ar>>alps::make_pvp("/"+cluster_name_+"/cluster_momenta", cluster_momenta_);
    ar>>alps::make_pvp("/"+cluster_name_+"/lattice_momenta", lattice_momenta_);
    ar>>alps::make_pvp("/"+cluster_name_+"/reciprocal_basis", reciprocal_basis_);
    ar>>alps::make_pvp("/"+cluster_name_+"/reciprocal_basis_inv", reciprocal_basis_inv_);
    ar>>alps::make_pvp("/"+cluster_name_+"/rbzone_boundaries", rbzone_boundaries_);
    ar>>alps::make_pvp("/"+cluster_name_+"/symmetry_table_k_space", symmetry_table_k_space_);
    ar>>alps::make_pvp("/"+cluster_name_+"/symmetry_table_real_space", symmetry_table_real_space_);
    
    read_site_class_set(ar, site_class_set_real_space_, "site_class_set_real_space",cluster_name_);
    read_site_class_set(ar, site_class_set_k_space_, "site_class_set_k_space",cluster_name_);
    read_pair_class_set(ar, pair_class_set_real_space_, "pair_class_set_real_space",cluster_name_);
    read_pair_class_set(ar, pair_class_set_k_space_, "pair_class_set_k_space",cluster_name_);
    return true; //success
}

std::vector<ClusterHelper::Matrix<double>::Type > ClusterHelper::setup_point_groups() const{
  //setup point group. This function has the 1d, 2d, and 3d simple cubic point groups hard coded.
  const int gr1D[] = { 1, -1 };
  const int gr2D[] = {  1,  0,  0,  1,
    -1,  0,  0,  1,
    1,  0,  0, -1,
    -1,  0,  0, -1,
    0,  1,  1,  0,
    0, -1,  1,  0,
    0,  1,  -1,  0,
    0, -1,  -1,  0 };
  const int gr3D[] = { 1, 0, 0,  0, 1, 0,  0, 0, 1,
    1, 0, 0,  0, -1, 0,  0, 0, -1,
    -1, 0, 0,  0, 1, 0,  0, 0, -1,
    -1, 0, 0,  0, -1, 0,  0, 0, -1,
    0, 0, 1,  1, 0, 0,  0, 1, 0,
    0, 0, -1,  -1, 0, 0,  0, 1, 0,
    0, 0, 1,  -1, 0, 0,  0, -1, 0,
    0, 0, -1,  1, 0, 0,  0, -1, 0,
    0, 1, 0,  0, 0, 1,  1, 0, 0,
    0, -1, 0,  0, 0, 1,  -1, 0, 0,
    0, -1, 0,  0, 0, -1,  1, 0, 0,
    0, 1, 0,  0, 0, -1,  -1, 0, 0,
    1, 0, 0,  0, 0, -1,  0, 1, 0,
    1, 0, 0,  0, 0, 1,  0, -1, 0,
    -1, 0, 0,  0, 0, -1,  0, -1, 0,
    -1, 0, 0,  0, 0, 1,  0, 1, 0,
    0, 0, 1,  0, 1, 0,  -1, 0, 0,
    0, 0, -1,  0, 1, 0,  1, 0, 0,
    0, 0, -1,  0, -1, 0,  -1, 0, 0,
    0, 0, 1,  0, -1, 0,  1, 0, 0,
    0, -1, 0,  1, 0, 0,  0, 0, 1,
    0, 1, 0,  -1, 0, 0,  0, 0, 1,
    0, -1, 0,  -1, 0, 0,  0, 0, -1,
    0, 1, 0,  1, 0, 0,  0, 0, -1,
    -1, 0, 0,  0, -1, 0,  0, 0, -1,
    -1, 0, 0,  0, 1, 0,  0, 0, 1,
    1, 0, 0,  0, -1, 0,  0, 0, 1,
    1, 0, 0,  0, 1, 0,  0, 0, -1,
    0, 0, -1,  -1, 0, 0,  0, -1, 0,
    0, 0, 1,  1, 0, 0,  0, -1, 0,
    0, 0, -1,  1, 0, 0,  0, 1, 0,
    0, 0, 1,  -1, 0, 0,  0, 1, 0,
    0, -1, 0,  0, 0, -1,  -1, 0, 0,
    0, 1, 0,  0, 0, -1,  1, 0, 0,
    0, 1, 0,  0, 0, 1,  -1, 0, 0,
    0, -1, 0,  0, 0, 1,  1, 0, 0,
    -1, 0, 0,  0, 0, 1,  0, -1, 0,
    -1, 0, 0,  0, 0, -1,  0, 1, 0,
    1, 0, 0,  0, 0, 1,  0, 1, 0,
    1, 0, 0,  0, 0, -1,  0, -1, 0,
    0, 0, -1,  0, -1, 0,  1, 0, 0,
    0, 0, 1,  0, -1, 0,  -1, 0, 0,
    0, 0, 1,  0, 1, 0,  1, 0, 0,
    0, 0, -1,  0, 1, 0,  -1, 0, 0,
    0, 1, 0,  -1, 0, 0,  0, 0, -1,
    0, -1, 0,  1, 0, 0,  0, 0, -1,
    0, 1, 0,  1, 0, 0,  0, 0, 1,
    0, -1, 0,  -1, 0, 0,  0, 0, 1 };
  const int* gr;
  unsigned int grsize;
  switch(dimension()) {
    case 1:
      gr = gr1D;
      grsize = 2;
      break;
    case 2:
      gr = gr2D;
      grsize = 8;
      break;
    case 3:
      gr = gr3D;
      grsize = 48;
      break;
    default:
      boost::throw_exception(std::logic_error("Only 1,2, or 3 dimensional clusters can be handled."));
      gr = gr1D;
      grsize = 0;
  }
  std::vector<Matrix<double>::Type > groups(grsize, Matrix<double>::Type(dimension(), dimension()));
  for (unsigned int g=0; g<grsize; ++g)
    for (unsigned int i=0; i<dimension(); ++i)
      for (unsigned int j=0; j<dimension(); ++j)
        groups[g](j,i) = gr[i+j*dimension()+g*dimension()*dimension()];
  return groups;
}
void ClusterHelper::symmetry_analysis(const Matrix<double>::Type& basis,
                                      const Matrix<double>::Type& basis_inv)
{
  //generate a table of all the point groups in this dimension
  std::vector<Matrix<double>::Type > groups=setup_point_groups();
  //setup symmetry table in real space
  for (unsigned int g=0; g<groups.size(); ++g) {
    std::map<int, int> vertexmap;
    std::set<int> vertexset;
    for (unsigned int i=0; i<n_cluster_sites(); ++i) {
      bool bad = false;
      for(int orig=0;orig<cell_origins().size(); ++orig){
        std::vector<double> c = cell_origins()[orig];
        for (unsigned int d=0; d<dimension(); ++d)
          c[d] += cluster_coordinate(i,d);
        c = matrix_vector_mult(groups[g], c);
        c = map_back_into_cluster(c, basis_inv, basis);
        unsigned int i_mapped = identify_vertex(c, cluster_coordinates_);
        if (!vertexmap.insert(std::pair<int,int>(i, i_mapped)).second) {
          if (vertexmap[i]!=i_mapped) {
            bad = true;
            break;
          }
        }
        else
          vertexset.insert(i_mapped);
      }
      if (bad)
        break;
    }
    if (vertexset.size()==n_cluster_sites() && vertexmap.size()==n_cluster_sites()) {
      append_column_no_duplicates(symmetry_table_real_space_, vertexmap);
    }
  }
  //setup equivalence classes of sites and site pairs
  boost::tie(site_class_set_real_space_, pair_class_set_real_space_) = setup_classes(symmetry_table_real_space_);
  //setup symmetry table in k space
  for (unsigned int g=0; g<groups.size(); ++g) {
    std::map<int, int> vertexmap;
    std::set<int> vertexset;
    for (unsigned int i=0; i<n_cluster_sites(); ++i) {
      bool bad = false;
      for(int off=0;off<cell_offsets().size();++off){
        std::vector<int> o = cell_offsets()[off];
        //std::cout << c[0] << " "<< c[1] << " mu ";
        //c = matrix_vector_mult(reciprocal_basis, c);
        //std::cout << o[0] << " "<< o[1] << " cm ";
        std::vector<double> c(dimension());
        for (unsigned int d=0; d<dimension(); ++d)
          c[d] = 2*M_PI*o[d]+cluster_momentum(i,d);
        //std::cout << c[0] << " "<< c[1] << " g ";
        c = matrix_vector_mult(groups[g], c);
        //std::cout << c[0] << " " << c[1] << " ma ";
        c = map_back_into_brillouin_zone(c);
        //std::cout << c[0] << " " << c[1] << " : ";
        unsigned int i_mapped = identify_vertex(c, cluster_momenta_);
        //std::cout << i << ":" << i_mapped << "\n";
        if (i_mapped == n_cluster_sites()) {
          bad = true;
          break;
        }
        if (!vertexmap.insert(std::pair<int,int>(i, i_mapped)).second) {
          if (vertexmap[i]!=i_mapped) {
            bad = true;
            break;
          }
        }
        else
          vertexset.insert(i_mapped);
      }
      if (bad)
        break;
    }
    if (vertexset.size()==n_cluster_sites() && vertexmap.size()==n_cluster_sites()) {
      /*if(alps::is_master())
       std::cout  << g << " ";*/
      append_column_no_duplicates(symmetry_table_k_space_, vertexmap);
    }
  }
  /*if(alps::is_master())
   std::cout << "\n";*/
  /*std::cout << "symmetry table k space:\n";
   for (unsigned int i=0; i<symmetry_table_k_space_.size1(); ++i) {
   for (unsigned int j=0; j<symmetry_table_k_space_.size2(); ++j) {
   std::cout << symmetry_table_k_space_(i,j) << " ";
   }
   std::cout << std::endl;
   }
   std::cout << "symmetry table k space (with momenta):\n";
   for (unsigned int i=0; i<symmetry_table_k_space_.size1(); ++i) {
   for (unsigned int j=0; j<symmetry_table_k_space_.size2(); ++j) {
   std::cout << " (";
   for(std::size_t d=0;d<dimension();++d)
   std::cout<<cluster_momentum(symmetry_table_k_space_(i,j),d) << " ";
   std::cout<<")";
   }
   std::cout << std::endl;
   }*/
  
  //setup equivalence classes of sites and site pairs
  boost::tie(site_class_set_k_space_, pair_class_set_k_space_) = setup_classes(symmetry_table_k_space_);
}



std::pair<ClusterHelper::site_class_set_type, ClusterHelper::pair_class_set_type>
ClusterHelper::setup_classes(const Matrix<int>::Type& symmetry_table) const
{
  //setup site equivalence classes
  site_class_set_type site_class_set;
  for (unsigned int j=0; j<symmetry_table.size2(); ++j) {
    site_class_type s_class;
    for (unsigned int i=0; i<symmetry_table.size1(); ++i) {
      s_class.insert(symmetry_table(i,j));
    }
    site_class_set.insert(s_class);
  }
  //setup site pair equivalence classes
  pair_class_set_type pair_class_set;
  for (unsigned int j1=0; j1<symmetry_table.size2(); ++j1) {
    for (unsigned int j2=0; j2<symmetry_table.size2(); ++j2) {
      pair_class_type p_class;
      for (unsigned int i=0; i<symmetry_table.size1(); ++i) {
        p_class.insert(pair_type(symmetry_table(i,j1), symmetry_table(i,j2)));
      }
      pair_class_set.insert(p_class);
    }
  }
  //reduce pair classes
  bool reducible = true;
  while (reducible) {
    reducible = false;
    pair_class_set_type inserted_class_set;
    pair_class_set_type reduced_pair_class_set;
    for (pair_class_set_type::const_iterator it_i1=pair_class_set.begin();
         it_i1!=pair_class_set.end(); ++it_i1) {
      if (inserted_class_set.find(*it_i1)==inserted_class_set.end()) {
        pair_class_type red_p_class = *it_i1;
        pair_class_set_type::const_iterator it_next = it_i1;
        ++it_next;
        for (pair_class_set_type::const_iterator it_i2=it_next;
             it_i2!=pair_class_set.end(); ++it_i2) {
          bool intersec = false;
          for (pair_class_type::const_iterator it_j1=it_i1->begin();
               it_j1!=it_i1->end(); ++it_j1) {
            for (pair_class_type::const_iterator it_j2=it_i2->begin();it_j2!=it_i2->end(); ++it_j2) {
              if (*it_j1==*it_j2 || (it_j1->first==it_j2->second && it_j1->second==it_j2->first&&paramagnet_)) {
                intersec = true;
                reducible = true;
                break;
              }
            }
            if (intersec)
              break;
          }
          if (intersec) {
            for (pair_class_type::const_iterator it_j2=it_i2->begin();
                 it_j2!=it_i2->end(); ++it_j2)
              red_p_class.insert(*it_j2);
            inserted_class_set.insert(*it_i2);
          }
        }
        reduced_pair_class_set.insert(red_p_class);
      }
    }
    pair_class_set = reduced_pair_class_set;
  }
  return std::pair<site_class_set_type, pair_class_set_type>(site_class_set, pair_class_set);
}


double ClusterHelper::site_distance(const unsigned int i, const unsigned int j) const
{
  double dist = 0;
  for (std::size_t d=0; d<dimension(); ++d) {
    double del = cluster_coordinate(i, d);
    del -= cluster_coordinate(j, d);
    dist += del*del;
  }
  return round(sqrt(dist));
}



bool ClusterHelper::inside_cluster(std::vector<double> v, const Matrix<double>::Type& B) const
{
  v = matrix_vector_mult(B, v);
  v = round(v);
  bool inside = true;
  for (unsigned int d=0; d<v.size(); ++d)
    inside &= v[d]<1 && v[d]>=0;
  return inside;
}



bool ClusterHelper::inside_reduced_brillouin_zone(std::vector<double> v) const
{
  for (std::size_t i=0; i<rbzone_boundaries_.size1(); ++i) {
    double norm = 0.;
    for (std::size_t d=0; d<rbzone_boundaries_.size2(); ++d)
      norm += rbzone_boundaries_(i,d)*rbzone_boundaries_(i,d);
    //norm = sqrt(norm);
    double projection = 0.;
    for (std::size_t d=0; d<rbzone_boundaries_.size2(); ++d)
      projection +=  v[d] * rbzone_boundaries_(i,d);
    projection /= norm;
    projection = round(projection);
    if (projection > 1 || projection <= -1)
      return false;
  }
  return true;
}



bool ClusterHelper::append_column_no_duplicates(Matrix<double>::Type& m, const std::vector<double>& v) const
{
  for (unsigned int i=0; i<m.size1(); ++i) {
    bool equal = true;
    for (unsigned int d=0; d<m.size2(); ++d)
      equal &= almost_equal(m(i,d), v[d]);
    if (equal)
      return false;
  }
  m.resize_preserve_values(m.size1()+1, v.size());
  for (unsigned int d=0; d<m.size2(); ++d)
    m(m.size1()-1,d) = v[d];
  return true;
}



bool ClusterHelper::append_column_no_duplicates(Matrix<int>::Type& m,
                                                std::map<int, int>& smap) const
{
  for (unsigned int i=0; i<m.size1(); ++i) {
    bool equal = true;
    for (unsigned int d=0; d<m.size2(); ++d)
      equal &= m(i,d)==smap[d];
    if (equal)
      return false;
  }
  m.resize_preserve_values(m.size1()+1, smap.size());
  for (unsigned int d=0; d<m.size2(); ++d)
    m(m.size1()-1,d) = smap[d];
  return true;
}



std::vector<double> ClusterHelper::map_back_into_brillouin_zone(std::vector<double> v) const
{
  for (unsigned int d=0; d<v.size(); ++d) {
    v[d] /= M_PI;
    v[d] = round(v[d]);
    while (v[d]>1) {
      v[d] -= 2;
    }
    while (v[d]<=-1) {
      v[d] += 2;
    }
    v[d] *= M_PI;
  }
  /*for (unsigned int d=0; d<v.size(); ++d) {
   std::cout << v[d];
   }std::cout<<std::endl;*/
  return v;
}

void ClusterHelper::write_hdf5(alps::hdf5::archive &ar, const std::string &p) const{
  ar<<alps::make_pvp(p+"/cluster_coordinates", cluster_coordinates_);
  ar<<alps::make_pvp(p+"/cluster_momenta", cluster_momenta_);
  ar<<alps::make_pvp(p+"/lattice_momenta", lattice_momenta_);
}


std::vector<double> ClusterHelper::map_back_into_cluster(std::vector<double> v,
                                                         const Matrix<double>::Type& B_inv,
                                                         const Matrix<double>::Type& B) const
{
  v = matrix_vector_mult(B_inv, v);
  for (unsigned int d=0; d<v.size(); ++d) {
    v[d] = round(v[d]);
    while (v[d]>=1) {
      v[d] -= 1;
    }
    while (v[d]<0) {
      v[d] += 1;
    }
  }
  v = matrix_vector_mult(B, v);
  return v;
}



int ClusterHelper::cluster_momentum_approximation(std::vector<double> v) const
{
  Matrix<double>::Type T;
  switch (dimension()) {
    case 1 :
      T.resize(1,1);
      T(0,0) = 2*M_PI;
      break;
    case 2 :
      T.resize(4,2);
      T(0,0) = 2*M_PI; T(0,1) = 0;
      T(1,0) = 0 ; T(1,1) = 2*M_PI;
      T(2,0) = 2*M_PI; T(2,1) = 2*M_PI;
      T(3,0) = 2*M_PI; T(3,1) = -2*M_PI;
      break;
    case 3 :
      T.resize(13,3);
      T(0,0)  =  2*M_PI; T(0,1)  =  0; T(0,2)  =  0;
      T(1,0)  =  0; T(1,1)  =  2*M_PI; T(1,2)  =  0;
      T(2,0)  =  0; T(2,1)  =  0; T(2,2)  =  2*M_PI;
      T(3,0)  =  2*M_PI; T(3,1)  =  2*M_PI; T(3,2)  =  0;
      T(4,0)  =  2*M_PI; T(4,1)  =  0; T(4,2)  =  2*M_PI;
      T(5,0)  =  0; T(5,1)  =  2*M_PI; T(5,2)  =  2*M_PI;
      T(6,0)  =  2*M_PI; T(6,1)  = -2*M_PI; T(6,2)  =  0;
      T(7,0)  =  2*M_PI; T(7,1)  =  0; T(7,2)  = -2*M_PI;
      T(8,0)  =  0; T(8,1)  =  2*M_PI; T(8,2)  = -2*M_PI;
      T(9,0)  =  2*M_PI; T(9,1)  =  2*M_PI; T(9,2)  =  2*M_PI;
      T(10,0) = -2*M_PI; T(10,1) =  2*M_PI; T(10,2) =  2*M_PI;
      T(11,0) =  2*M_PI; T(11,1) = -2*M_PI; T(11,2) =  2*M_PI;
      T(12,0) =  2*M_PI; T(12,1) =  2*M_PI; T(12,2) = -2*M_PI;
  }
  for (std::size_t i=0; i<cluster_momenta_.size1(); ++i) {
    for (std::size_t t=0; t<T.size1(); ++t) {
      for (int s=-1; s<=1; ++s) {
        std::vector<double> diff(dimension());
        for (unsigned int d=0; d<dimension(); ++d)
          diff[d] = v[d] - cluster_momenta_(i,d) - s*T(t,d);
        diff = matrix_vector_mult(reciprocal_basis_inv_, diff);
        for (unsigned int d=0; d<dimension(); ++d)
          diff[d] *= 2;
        diff = round(diff);
        bool inside = true;
        for (unsigned int d=0; d<diff.size(); ++d)
          inside &= diff[d]<=1 && diff[d]>-1;
        if (inside) {
          return i;
        }
      }
    }
  }
  boost::throw_exception(std::logic_error("no cluster momentum approximation found"));
  return n_cluster_sites();
}



unsigned int ClusterHelper::identify_vertex(const std::vector<double>& v, const Matrix<double>::Type& coordinate_matrix) const
{
  for (unsigned int i=0; i<coordinate_matrix.size1(); ++i) {
    bool equal = true;
    for (unsigned int d=0; d<coordinate_matrix.size2(); ++d)
      equal &= almost_equal(v[d], coordinate_matrix(i,d));
    if (equal)
      return i;
  }
  //std::cout << v[0] << " " << v[1] << "\n";
  return n_cluster_sites();
}


unsigned int ClusterHelper::identify_sublattice(const Matrix<double>::Type& coordinate_matrix, const unsigned int i) const
{
  double sum=0;
  for (std::size_t d=0; d<coordinate_matrix.size2(); ++d)
    sum += coordinate_matrix(i,d);
  sum = round(sum);
  sum *= sum<0 ? -1 : 1;
  return (int(round(sum)))%2;
}

double ClusterHelper::round(const double d) const
{
  int rounded = int(d > 0.0 ? d + 0.5 : d - 0.5);
  if (almost_equal(rounded, d))
    return rounded;
  else
    return d;
}



std::vector<double> ClusterHelper::round(std::vector<double> v) const
{
  for (unsigned int d=0; d<v.size(); ++d)
    v[d] = round(v[d]);
  return v;
}



std::vector<double> ClusterHelper::matrix_vector_mult(const Matrix<double>::Type& M, const std::vector<double>& v) const
{
  std::vector<double> vn(v.size());
  for (unsigned int i=0; i<v.size(); ++i) {
    vn[i] = 0.;
    for (unsigned int j=0; j<v.size(); ++j) {
      vn[i] += M(j,i)*v[j];
    }
  }
  return vn;
}

bool ClusterHelper::operator==(const ClusterHelper &CH) const{
  //compare scalars
  if(epsilon_!=CH.epsilon_) return false;
  if(bipartite_!=CH.bipartite_) return false;
  if(paramagnet_!=CH.paramagnet_) return false;
  
  //compare double matrices
  if(Matrix<double>::norm_1(cluster_coordinates_ - CH.cluster_coordinates_) > epsilon_) return false;
  if(Matrix<double>::norm_1(lattice_momenta_ - CH.lattice_momenta_) > epsilon_) return false;
  if(Matrix<double>::norm_1(reciprocal_basis_ - CH.reciprocal_basis_) > epsilon_) return false;
  if(Matrix<double>::norm_1(reciprocal_basis_inv_ - CH.reciprocal_basis_inv_) > epsilon_) return false;
  if(Matrix<double>::norm_1(rbzone_boundaries_ - CH.rbzone_boundaries_) > epsilon_) return false;
  
  //compare int matrices
  if(Matrix<int>::norm_1(symmetry_table_k_space_ - CH.symmetry_table_k_space_) > epsilon_) return false;
  if(Matrix<int>::norm_1(symmetry_table_real_space_ - CH.symmetry_table_real_space_) > epsilon_) return false;
  
  for(int i=0;i<periodicity_.size();++i){
    if(periodicity_[i]!=CH.periodicity_[i]) return false;
  }
  
  return true;
  
}
