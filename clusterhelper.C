/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2007 by Sebastian Fuchs <fuchs@comp-phys.org>
 * Copyright (C) 2007-2015 by Emanuel Gull <egull@umich.edu>
 *
 *
 *****************************************************************************/

#include "alpsgraph.h"
#include "cluster.h"
#include "clusterhelper.hxx"
#include "alps/hdf5/ublas/vector.hpp"
#include "alps/hdf5/ublas/matrix.hpp"


ClusterHelper::ClusterHelper(const alps::params& params) :
alps_graph(params["LATTICE"],params["LATTICE_LIBRARY"],params["L"]),
paramagnet_(params["PARAMAGNET"]),
epsilon_(1e-9)
{

  if (int(params["SITES"]) != alps_graph::unit_cell_volume()) throw std::runtime_error("the graph"+params["LATTICE"]+" unit cell volume is not :"+params["SITES"]);
  //try to read this transformer, otherwise generate it.
  if(!read(params))
    setup_tables(params);

  //Find the periodicity in each direction
  //Set Nxyz
  	Nxyz.clear();
  	Nxyz.resize(dimension(),0.);
  	for(int k=0; k<dimension(); k++){
  		int n=0;
  		bool done=false;
  		while(!done){
  			n+=1;
  			for(int d=0; d<dimension(); d++){
  				if(n*reciprocal_basis_(k,d) > 2*M_PI){
  					done = true;
  				}
  			}
  		}
  		Nxyz[k] = n-1;
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

void ClusterHelper::setup_tables(const alps::params &params)
{
  //std::cout<<"setting up tables of transformer"<<std::endl;
  //setup basis vectors
  Matrix<double>::Type basis(dimension(), dimension());
  memset(&(basis(0,0)), 0, sizeof(double)*dimension()*dimension());
  //dca_graph::basis_vector_iterator itb, itb_end;
  //boost::tie(itb, itb_end) = basis_vector_iterators();
  int sum = 0;
  for (unsigned int i=0; i<basis_vectors().size(); ++i) {
    int sum_i = 0;
    for (unsigned int d=0; d<dimension(); ++d) {
      basis(i,d) = basis_vectors()[i][d];
      sum_i += int(round(basis(i,d)));
    }
    sum += sum_i>=0 ? sum_i%2 : (-sum_i)%2;
  }
  bipartite_ = (sum==0);
  Matrix<double>::Type basis_inv = Matrix<double>::inverse(basis);
  //setup cluster coordinates
  int maxc=n_cluster_sites();
  std::vector<double> c(dimension());
  switch(dimension()) {
    case 1 :
      for (int i=-maxc; i<=maxc; ++i) {
        c[0] = i;
        if (inside_cluster(c, basis_inv))
          append_column_no_duplicates(cluster_coordinates_, c);
      }
      break;
    case 2:
      for (int i=-maxc; i<=maxc; ++i) {
        for (int j=-maxc; j<=maxc; ++j) {
          c[0] = i; c[1] = j;
          if (inside_cluster(c, basis_inv))
            append_column_no_duplicates(cluster_coordinates_, c);
        }
      }
      break;
    case 3:
      for (int i=-maxc; i<=maxc; ++i) {
        for (int j=-maxc; j<=maxc; ++j) {
          for (int k=-maxc; k<=maxc; ++k) {
            c[0] = i; c[1] = j; c[2]=k;
            if (inside_cluster(c, basis_inv)) {
              append_column_no_duplicates(cluster_coordinates_, c);
            }
          }
        }
      }
      break;
  }
  //arrange vertices for bipartite lattices
  if (bipartite()) {
    /*if(alps::is_master())
     std::cout << "cluster is bipartite\n";*/
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
  /*if(alps::is_master()){
    std::cerr<<"DEBUG: Writing cluster coordinates"<<std::endl;
    std::ofstream clustercoordstream("ClusterCoordinates");
    for (std::size_t m=0; m<cluster_coordinates_.size1(); ++m) {
      clustercoordstream << m << " ";
      for (std::size_t d=0; d<cluster_coordinates_.size2(); ++d)
        clustercoordstream << cluster_coordinates_(m,d) << " ";
      clustercoordstream << std::endl;
    }
  }*/
  //setup reciprocal basis
  reciprocal_basis_.resize(dimension(), dimension());
  memset(&(reciprocal_basis_(0,0)), 0, sizeof(double)*dimension()*dimension());
  for (unsigned int i=0; i<reciprocal_basis_vectors().size();++i)
    for (unsigned int d=0; d<dimension(); ++d)
      reciprocal_basis_(i,d) = reciprocal_basis_vectors()[i][d];
  reciprocal_basis_inv_ = Matrix<double>::inverse(reciprocal_basis_);
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
  symmetry_analysis(basis, basis_inv/*, reciprocal_basis, reciprocal_basis_inv*/);
#ifndef USE_MPI
  write(params);
#else
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0){
    write(params);
  }
#endif
}
void ClusterHelper::write_site_class_set(alps::hdf5::archive &ar, const site_class_set_type &set, const std::string &name) const{
  int z=0;
  std::stringstream path; path<<"/transformer/"<<name<<"/size";
  ar<<alps::make_pvp(path.str(),set.size());
  for(site_class_set_type::const_iterator it=set.begin(); it != set.end();++it){
    int z1=0;
    std::stringstream path; path<<"/transformer/"<<name<<"/"<<z<<"/size";
    ar<<alps::make_pvp(path.str(),it->size());
    for(site_class_type::const_iterator it1=it->begin(); it1 != it->end();++it1){
      std::stringstream path; path<<"/transformer/"<<name<<"/"<<z<<"/"<<z1;
      ar<<alps::make_pvp(path.str(), *it1);
      z1++;
    }
    z++;
  }
}
void ClusterHelper::write_pair_class_set(alps::hdf5::archive &ar, const pair_class_set_type &set, const std::string &name) const{
  int z=0;
  std::stringstream path; path<<"/transformer/"<<name<<"/size";
  ar<<alps::make_pvp(path.str(),set.size());
  for(pair_class_set_type::const_iterator it=set.begin(); it != set.end();++it){
    int z1=0;
    std::stringstream path; path<<"/transformer/"<<name<<"/"<<z<<"/size";
    ar<<alps::make_pvp(path.str(),it->size());
    for(pair_class_type::const_iterator it1=it->begin(); it1 != it->end();++it1){
      std::stringstream path1, path2;
      path1<<"/transformer/"<<name<<"/"<<z<<"/"<<z1<<"/first";
      ar<<alps::make_pvp(path1.str(), it1->first);
      path2<<"/transformer/"<<name<<"/"<<z<<"/"<<z1<<"/second";
      ar<<alps::make_pvp(path2.str(), it1->second);
      z1++;
    }
    z++;
  }
}
void ClusterHelper::read_site_class_set(alps::hdf5::archive &ar, site_class_set_type &set, const std::string &name){
  int zmax, z1max;
  std::stringstream path; path<<"/transformer/"<<name<<"/size";
  ar>>alps::make_pvp(path.str(),zmax);
  for(int z=0;z<zmax;++z){
    std::stringstream path; path<<"/transformer/"<<name<<"/"<<z<<"/size";
    ar>>alps::make_pvp(path.str(),z1max);
    std::set<site_t> s;
    for(int z1=0;z1<z1max;++z1){
      site_t val;
      std::stringstream path; path<<"/transformer/"<<name<<"/"<<z<<"/"<<z1;
      ar>>alps::make_pvp(path.str(), val);
      s.insert(val);
    }
    set.insert(s);
  }
}
void ClusterHelper::read_pair_class_set(alps::hdf5::archive &ar, pair_class_set_type &set, const std::string &name){
  int zmax, z1max;
  std::stringstream path; path<<"/transformer/"<<name<<"/size";
  ar>>alps::make_pvp(path.str(),zmax);
  for(int z=0;z<zmax;++z){
    std::stringstream path; path<<"/transformer/"<<name<<"/"<<z<<"/size";
    ar>>alps::make_pvp(path.str(),z1max);
    std::set<pair_type> s;
    for(int z1=0;z1<z1max;++z1){
      pair_type val;
      std::stringstream path1; path1<<"/transformer/"<<name<<"/"<<z<<"/"<<z1<<"/first";
      std::stringstream path2; path2<<"/transformer/"<<name<<"/"<<z<<"/"<<z1<<"/second";
      ar>>alps::make_pvp(path1.str(), val.first);
      ar>>alps::make_pvp(path2.str(), val.second);
      s.insert(val);
    }
    set.insert(s);
  }
}

bool ClusterHelper::write(const alps::params &parms) const{
  ///TODO: this should implement writing the clusterhelper.
  if(parms.defined("TMPNAME")){
    //write data to the 'tmpname' file
    alps::hdf5::archive ar(parms["TMPNAME"]+std::string(".in.h5"), alps::hdf5::archive::WRITE);
    if(ar.is_data("/transformer/name")) return false; //already have a transformer there.
    std::string latticename=parms["LATTICE"];
    ar<<alps::make_pvp("/transformer/name", latticename);
    ar<<alps::make_pvp("/transformer/bipartite", bipartite_);
    ar<<alps::make_pvp("/transformer/epsilon", epsilon_);
    ar<<alps::make_pvp("/transformer/cluster_coordinates", cluster_coordinates_);
    ar<<alps::make_pvp("/transformer/cluster_momenta", cluster_momenta_);
    ar<<alps::make_pvp("/transformer/lattice_momenta", lattice_momenta_);
    ar<<alps::make_pvp("/transformer/reciprocal_basis", reciprocal_basis_);
    ar<<alps::make_pvp("/transformer/reciprocal_basis_inv", reciprocal_basis_inv_);
    ar<<alps::make_pvp("/transformer/rbzone_boundaries", rbzone_boundaries_);
    ar<<alps::make_pvp("/transformer/symmetry_table_k_space", symmetry_table_k_space_);
    ar<<alps::make_pvp("/transformer/symmetry_table_real_space", symmetry_table_real_space_);
    
    //the sets we need to do by hand:
    //site class set real space
    write_site_class_set(ar, site_class_set_real_space_, "site_class_set_real_space");
    write_site_class_set(ar, site_class_set_k_space_, "site_class_set_k_space");
    write_pair_class_set(ar, pair_class_set_real_space_, "pair_class_set_real_space");
    write_pair_class_set(ar, pair_class_set_k_space_, "pair_class_set_k_space");
    return true;
  }
  return false;
}
bool ClusterHelper::read(const alps::params &parms){
  ///TODO: this should implement writing the clusterhelper.
  if(parms.defined("TMPNAME")){
    {
      std::string tmpname=parms["TMPNAME"];
      std::ifstream in_file(tmpname.c_str());
      if(!in_file.is_open()) return false;
    }
    alps::hdf5::archive ar(parms["TMPNAME"]+std::string(".in.h5"), alps::hdf5::archive::READ);
    if(!ar.is_data("/transformer/name")) return false;
    std::string name;     ar>>alps::make_pvp("/transformer/name", name);
    std::string latticename=parms["LATTICE"];
    if(name!=latticename) throw std::logic_error("reading cluster transformer with name: "+name+" for cluster: "+latticename);
    //write data to the 'tmpname' file
    ar>>alps::make_pvp("/transformer/bipartite", bipartite_);
    ar>>alps::make_pvp("/transformer/epsilon", *(const_cast<double*>(&epsilon_)));
    ar>>alps::make_pvp("/transformer/cluster_coordinates", cluster_coordinates_);
    ar>>alps::make_pvp("/transformer/cluster_momenta", cluster_momenta_);
    ar>>alps::make_pvp("/transformer/lattice_momenta", lattice_momenta_);
    ar>>alps::make_pvp("/transformer/reciprocal_basis", reciprocal_basis_);
    ar>>alps::make_pvp("/transformer/reciprocal_basis_inv", reciprocal_basis_inv_);
    ar>>alps::make_pvp("/transformer/rbzone_boundaries", rbzone_boundaries_);
    ar>>alps::make_pvp("/transformer/symmetry_table_k_space", symmetry_table_k_space_);
    ar>>alps::make_pvp("/transformer/symmetry_table_real_space", symmetry_table_real_space_);
    
    read_site_class_set(ar, site_class_set_real_space_, "site_class_set_real_space");
    read_site_class_set(ar, site_class_set_k_space_, "site_class_set_k_space");
    read_pair_class_set(ar, pair_class_set_real_space_, "pair_class_set_real_space");
    read_pair_class_set(ar, pair_class_set_k_space_, "pair_class_set_k_space");
    return true; //success
  }
  return false;
}

void ClusterHelper::symmetry_analysis(const Matrix<double>::Type& basis,
                                      const Matrix<double>::Type& basis_inv/*,
                                                                            const Matrix<double>::Type& reciprocal_basis,
                                                                            const Matrix<double>::Type& reciprocal_basis_inv */)
{
  //setup point group
  int gr1D[] = { 1, -1 };
  int gr2D[] = {  1,  0,  0,  1,
    -1,  0,  0,  1,
    1,  0,  0, -1,
    -1,  0,  0, -1,
    0,  1,  1,  0,
    0, -1,  1,  0,
    0,  1,  -1,  0,
    0, -1,  -1,  0 };
  int gr3D[] = { 1, 0, 0,  0, 1, 0,  0, 0, 1,
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
  int* gr;
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
  std::vector<Matrix<double>::Type > group(grsize, Matrix<double>::Type(dimension(), dimension()));
  for (unsigned int g=0; g<grsize; ++g)
    for (unsigned int i=0; i<dimension(); ++i)
      for (unsigned int j=0; j<dimension(); ++j)
        group[g](j,i) = gr[i+j*dimension()+g*dimension()*dimension()];
  //setup symmetry table in real space
  //Matrix<site_t>::Type symmetry_table_real_space;
  /*if(alps::is_master())
   std::cout << "good symmetries: ";*/
  for (unsigned int g=0; g<group.size(); ++g) {
    std::map<site_t, site_t> vertexmap;
    std::set<site_t> vertexset;
    for (unsigned int i=0; i<n_cluster_sites(); ++i) {
      bool bad = false;
      //-> at this point we should separate knowledge of cells from knowledge of the alps graph structure. All we need is an origin, which is a vector type.
      //we can achieve this by generating a vector of origins and iterating over them.
      for(int orig=0;orig<cell_origins().size(); ++orig){
        std::vector<double> c = cell_origins()[orig];
        for (unsigned int d=0; d<dimension(); ++d)
          c[d] += cluster_coordinate(i,d);
        c = matrix_vector_mult(group[g], c);
        c = map_back_into_cluster(c, basis_inv, basis);
        unsigned int i_mapped = identify_vertex(c, cluster_coordinates_);
        if (!vertexmap.insert(std::pair<site_t,site_t>(i, i_mapped)).second) {
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
      append_column_no_duplicates(symmetry_table_real_space_, vertexmap);
    }
  }
  //This part is for foursite only, needs generalization
  {
    std::map<site_t, site_t> vertexmap;
    
    vertexmap[0]=2;vertexmap[1]=1;vertexmap[2]=0;vertexmap[3]=3;
    append_column_no_duplicates(symmetry_table_real_space_, vertexmap);
    vertexmap[0]=2;vertexmap[1]=3;vertexmap[2]=0;vertexmap[3]=1;
    append_column_no_duplicates(symmetry_table_real_space_, vertexmap);
  }
 // std::cout << "\n";
 /* std::cout << "symmetry table in real space:\n";
   for (unsigned int i=0; i<symmetry_table_real_space_.size1(); ++i) {
   for (unsigned int j=0; j<symmetry_table_real_space_.size2(); ++j) {
   std::cout << symmetry_table_real_space_(i,j) << " ";
   }
   std::cout << std::endl;
   }*/
  //setup equivalence classes of sites and site pairs
  boost::tie(site_class_set_real_space_, pair_class_set_real_space_) = setup_classes(symmetry_table_real_space_);
  //setup symmetry table in k space
  //Matrix<site_t>::Type symmetry_table_k_space_;
  /*if(alps::is_master())
   std::cout << "good symmetries: ";*/
  for (unsigned int g=0; g<group.size(); ++g) {
    std::map<site_t, site_t> vertexmap;
    std::set<site_t> vertexset;
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
        c = matrix_vector_mult(group[g], c);
        //std::cout << c[0] << " " << c[1] << " ma ";
        c = map_back_into_brillouin_zone(c);
        //std::cout << c[0] << " " << c[1] << " : ";
        unsigned int i_mapped = identify_vertex(c, cluster_momenta_);
        //std::cout << i << ":" << i_mapped << "\n";
        if (i_mapped == n_cluster_sites()) {
          bad = true;
          break;
        }
        if (!vertexmap.insert(std::pair<site_t,site_t>(i, i_mapped)).second) {
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
ClusterHelper::setup_classes(const Matrix<site_t>::Type& symmetry_table) const
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
  m.resize(m.size1()+1, v.size(), true);
  for (unsigned int d=0; d<m.size2(); ++d)
    m(m.size1()-1,d) = v[d];
    return true;
}



bool ClusterHelper::append_column_no_duplicates(Matrix<site_t>::Type& m,
                                                std::map<site_t, site_t>& smap) const
{
  for (unsigned int i=0; i<m.size1(); ++i) {
    bool equal = true;
    for (unsigned int d=0; d<m.size2(); ++d)
      equal &= m(i,d)==smap[d];
    if (equal)
      return false;
  }
  m.resize(m.size1()+1, smap.size(), true);
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

