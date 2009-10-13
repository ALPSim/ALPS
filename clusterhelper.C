/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2007 by Sebastian Fuchs <fuchs@comp-phys.org>
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


#include "cluster.h"
#include<alps/osiris/comm.h>


ClusterHelper::ClusterHelper(const alps::Parameters& params) :
  alps::graph_helper<alps::graph_helper<>::graph_type>(params) ,  
  n_cluster_sites_(boost::lexical_cast<int>(params["SITES"])), epsilon_(1e-9)
{
  setup_tables();
  /*std::ofstream siterealstr("SiteClassesRealSpace");
  for (std::set<site_class_type>::const_iterator it_i=site_class_set_real_space_.begin(); 
       it_i!=site_class_set_real_space_.end(); ++it_i) {
    for (site_class_type::const_iterator it_j=it_i->begin(); it_j!=it_i->end(); ++it_j) {
      siterealstr << *it_j << " ";
    }
    siterealstr << std::endl;
  }
  std::ofstream pairrealstr("PairClassesRealSpace");
  for (std::set<pair_class_type >::const_iterator it_i=pair_class_set_real_space_.begin(); 
       it_i!=pair_class_set_real_space_.end(); ++it_i) {
    for (pair_class_type::const_iterator it_j=it_i->begin(); it_j!=it_i->end(); ++it_j) {
      pairrealstr << "(" << it_j->first << "," << it_j->second << ") ";
    }
    pairrealstr << std::endl;
  }
  std::ofstream sitekstr("SiteClassesKSpace");
  for (std::set<site_class_type>::const_iterator it_i=site_class_set_k_space_.begin(); 
       it_i!=site_class_set_k_space_.end(); ++it_i) {
    for (site_class_type::const_iterator it_j=it_i->begin(); it_j!=it_i->end(); ++it_j) {
      sitekstr << *it_j << " ";
    }
    sitekstr << std::endl;
  }
  std::ofstream pairkstr("PairClassesKSpace");
  for (std::set<pair_class_type >::const_iterator it_i=pair_class_set_k_space_.begin(); 
       it_i!=pair_class_set_k_space_.end(); ++it_i) {
    for (pair_class_type::const_iterator it_j=it_i->begin(); it_j!=it_i->end(); ++it_j) {
      pairkstr << "(" << it_j->first << "," << it_j->second << ") ";
    }
    pairkstr << std::endl;
    }*/
}



void ClusterHelper::setup_tables()
{ 
  //setup basis vectors
  Matrix<double>::Type basis(dimension(), dimension());
  basis_vector_iterator itb, itb_end;
  boost::tie(itb, itb_end) = basis_vectors();
  int sum = 0;
  for (unsigned int i=0; itb!=itb_end; ++itb, ++i) {
    int sum_i = 0;
    for (unsigned int d=0; d<dimension(); ++d) {
      basis(i,d) = (*itb)[d];
      sum_i += int(round(basis(i,d)));
    }
    sum += sum_i>=0 ? sum_i%2 : (-sum_i)%2;
  }
  bipartite_ = (sum==0);
  Matrix<double>::Type basis_inv = Matrix<double>::inverse(basis);
  //setup cluster coordinates
  int maxc=n_cluster_sites_;
  vector_type c(dimension());
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
    std::cout << "cluster is bipartite\n";
    for (int i=0; i<cluster_coordinates_.size1(); ++i) {
      if (identify_sublattice(cluster_coordinates_, i)!=i%2) {
	int j = i;
	while (identify_sublattice(cluster_coordinates_, ++j)!=i%2) {};
	for (int d=0; d<cluster_coordinates_.size2(); ++d) {
	  double t=cluster_coordinates_(i,d);
	  cluster_coordinates_(i,d) = cluster_coordinates_(j,d);
	  cluster_coordinates_(j,d) = t;
	} 
      }
    }
  }
  if(alps::is_master()){
    std::ofstream clustercoordstream("ClusterCoordinates");
    for (int m=0; m<cluster_coordinates_.size1(); ++m) {
      clustercoordstream << m << " ";
      for (int d=0; d<cluster_coordinates_.size2(); ++d) 
	clustercoordstream << cluster_coordinates_(m,d) << " ";
      clustercoordstream << std::endl;
    }
  }
  //setup reciprocal basis
  reciprocal_basis_.resize(dimension(), dimension());
  basis_vector_iterator itbr, itbr_end;
  boost::tie(itbr, itbr_end) = reciprocal_basis_vectors();
  for (unsigned int i=0; itbr!=itbr_end; ++itbr, ++i) 
    for (unsigned int d=0; d<dimension(); ++d)
      reciprocal_basis_(i,d) = (*itbr)[d];
  reciprocal_basis_inv_ = Matrix<double>::inverse(reciprocal_basis_);
  //setup cluster momenta
  if (bipartite()) {
    switch (dimension()) {
    case 1 :
      rbzone_boundaries_.resize(1,1);
      rbzone_boundaries_(0,0) = M_PI/2;
      break;
    case 2 :
      rbzone_boundaries_.resize(2,2);
      rbzone_boundaries_(0,0) = M_PI/2; rbzone_boundaries_(0,1) = M_PI/2;
      rbzone_boundaries_(1,0) = -M_PI/2; rbzone_boundaries_(1,1) = M_PI/2;
      break;
    case 3 :
      //reduced Brillouin zone is the one of an fcc lattice with a=2
      rbzone_boundaries_.resize(4,3);
      rbzone_boundaries_(0,0) =  M_PI/2; rbzone_boundaries_(0,1) =  M_PI/2; rbzone_boundaries_(0,2) =  M_PI/2;
      rbzone_boundaries_(1,0) = -M_PI/2; rbzone_boundaries_(1,1) =  M_PI/2; rbzone_boundaries_(1,2) =  M_PI/2;
      rbzone_boundaries_(2,0) =  M_PI/2; rbzone_boundaries_(2,1) = -M_PI/2; rbzone_boundaries_(2,2) =  M_PI/2;
      rbzone_boundaries_(3,0) =  M_PI/2; rbzone_boundaries_(3,1) =  M_PI/2; rbzone_boundaries_(3,2) = -M_PI/2;
      break;
    }
    cell_iterator itc, itc_end;
    for (boost::tie(itc, itc_end) = cells(); itc!=itc_end; ++itc) {
      vector_type o = origin(*itc);
      for (unsigned int i=0; i<n_cluster_sites(); ++i) {
	vector_type c = o;
	for (unsigned int d=0; d<dimension(); ++d) 
	  c[d] += cluster_coordinate(i,d);
	c = matrix_vector_mult(reciprocal_basis_, c);
	c = map_back_into_brillouin_zone(c);
	if (inside_reduced_brillouin_zone(c)) {
	  append_column_no_duplicates(cluster_momenta_, c);
	  for (int d=0; d<dimension(); ++d) 
	    c[d] += M_PI;
	  c = map_back_into_brillouin_zone(c);
	  append_column_no_duplicates(cluster_momenta_, c);
	}
      }
    }
  }
  else {
    cell_iterator itc, itc_end;
    for (boost::tie(itc, itc_end) = cells(); itc!=itc_end; ++itc) {
      vector_type o = origin(*itc);
      for (unsigned int i=0; i<n_cluster_sites(); ++i) {
	vector_type c = o;
	for (unsigned int d=0; d<dimension(); ++d) 
	  c[d] += cluster_coordinate(i,d);
	c = matrix_vector_mult(reciprocal_basis_, c);
	c = map_back_into_brillouin_zone(c);
	append_column_no_duplicates(cluster_momenta_, c);
      }
    }
    if ((n_cluster_sites()%2)==0) {
      Matrix<double>::Type sorted_momenta;
      vector_type v(dimension());
      unsigned int ind = 0;
      for(unsigned int i=0; i<n_cluster_sites(); ++i) {
	//std::cout << ind << "\n";
	for (unsigned int d=0; d<dimension(); ++d) 
	  v[d] = cluster_momenta_(ind,d);
	append_column_no_duplicates(sorted_momenta, v);
	for (unsigned int d=0; d<dimension(); ++d) 
	  v[d] += M_PI;
	v = map_back_into_brillouin_zone(v);
	ind = cluster_momentum_approximation(v);
      }
      cluster_momenta_ = sorted_momenta;
    }
  }
  if(alps::is_master()){
    std::ofstream clustermomstream("ClusterMomenta");
    for (int m=0; m<cluster_momenta_.size1(); ++m) {
      clustermomstream << m << " ";
      for (int d=0; d<cluster_momenta_.size2(); ++d) 
	clustermomstream << cluster_momenta_(m,d) << " ";
	clustermomstream << std::endl;
    }
    }
  //setup lattice momenta
  lattice_momenta_.resize(n_cluster_cells(), dimension());
  momentum_iterator itm, itm_end;
  boost::tie(itm, itm_end) = momenta();
  for (unsigned int m=0; itm!=itm_end; ++itm, ++m) {
    vector_type mom=momentum(*itm);
    for (int d=0; d<dimension(); ++d) {
      lattice_momenta_(m,d) = mom[d];
    }
  } 
  /*std::ofstream latticemomstream("LatticeMomenta");
  for (int m=0; m<lattice_momenta_.size1(); ++m) {
    latticemomstream << m << " ";
    for (int d=0; d<lattice_momenta_.size2(); ++d) 
      latticemomstream << lattice_momenta_(m,d) << " ";
    latticemomstream << std::endl;
    }*/
  if (cluster_coordinates_.size1() != n_cluster_sites())
    boost::throw_exception(std::logic_error("parameter SITES does not match cluster definition"));
  if (cluster_momenta_.size1() != n_cluster_sites())
    boost::throw_exception(std::logic_error("setup of cluster momenta table failed"));
  symmetry_analysis(basis, basis_inv/*, reciprocal_basis, reciprocal_basis_inv*/);
}





void ClusterHelper::symmetry_analysis(const Matrix<double>::Type& basis, 
				      const Matrix<double>::Type& basis_inv/*, 
				      const Matrix<double>::Type& reciprocal_basis, 
				      const Matrix<double>::Type& reciprocal_basis_inv */)
{
  //setup point group
  int gr1D[] = { 1, -1 };
  int gr2D[] = {  1,  0,    0,  1,
		  -1,  0,    0,  1,
		  1,  0,    0, -1,
		  -1,  0,    0, -1,
		  0,  1,    1,  0,
		  0, -1,    1,  0,
		  0,  1,   -1,  0,
		  0, -1,   -1,  0 };
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
  Matrix<site_t>::Type symmetry_table_real_space;
  std::cout << "good symmetries: ";
  for (unsigned int g=0; g<group.size(); ++g) {
    std::map<site_t, site_t> vertexmap;
    std::set<site_t> vertexset;
    for (unsigned int i=0; i<n_cluster_sites(); ++i) {
      bool bad = false;
      cell_iterator itc, itc_end;
      for (boost::tie(itc, itc_end) = cells(); itc!=itc_end; ++itc) {
	vector_type c = origin(*itc);
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
      std::cout  << g << " ";
      append_column_no_duplicates(symmetry_table_real_space, vertexmap);
    }
  }
  std::cout << "\n";
  /*std::cout << "symmetry table in real space:\n";
  for (unsigned int i=0; i<symmetry_table_real_space.size1(); ++i) {
    for (unsigned int j=0; j<symmetry_table_real_space.size2(); ++j) {
      std::cout << symmetry_table_real_space(i,j) << " ";
    }
    std::cout << std::endl;
    }*/
  //setup equivalence classes of sites and site pairs
  boost::tie(site_class_set_real_space_, pair_class_set_real_space_) = setup_classes(symmetry_table_real_space);
  //setup symmetry table in k space
  Matrix<site_t>::Type symmetry_table_k_space;
  std::cout << "good symmetries: ";
  for (unsigned int g=0; g<group.size(); ++g) {
    std::map<site_t, site_t> vertexmap;
    std::set<site_t> vertexset;
    for (unsigned int i=0; i<n_cluster_sites(); ++i) {
      bool bad = false;
      cell_iterator itc, itc_end;
      for (boost::tie(itc, itc_end) = cells(); itc!=itc_end; ++itc) {
	offset_type o = offset(*itc);
	//std::cout << c[0] << " "<< c[1] << " mu ";
	//c = matrix_vector_mult(reciprocal_basis, c);
	//std::cout << o[0] << " "<< o[1] << " cm ";
	vector_type c(dimension());
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
      std::cout  << g << " ";
      append_column_no_duplicates(symmetry_table_k_space, vertexmap);
    }
  }
  std::cout << "\n";
  /*std::cout << "symmetry table k space:\n";
  for (unsigned int i=0; i<symmetry_table_k_space.size1(); ++i) {
    for (unsigned int j=0; j<symmetry_table_k_space.size2(); ++j) {
      std::cout << symmetry_table_k_space(i,j) << " ";
    }
    std::cout << std::endl;
    }*/
  //setup equivalence classes of sites and site pairs
  boost::tie(site_class_set_k_space_, pair_class_set_k_space_) = setup_classes(symmetry_table_k_space);
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
	    for (pair_class_type::const_iterator it_j2=it_i2->begin(); 
		 it_j2!=it_i2->end(); ++it_j2) {
	      if (*it_j1==*it_j2 || (it_j1->first==it_j2->second && it_j1->second==it_j2->first)) {
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



bool ClusterHelper::inside_cluster(vector_type v, const Matrix<double>::Type& B) const 
{
  v = matrix_vector_mult(B, v);
  v = round(v);
  bool inside = true;
  for (unsigned int d=0; d<v.size(); ++d) 
    inside &= v[d]<1 && v[d]>=0;
  return inside;
}



bool ClusterHelper::inside_reduced_brillouin_zone(vector_type v) const 
{
  for (int i=0; i<rbzone_boundaries_.size1(); ++i) {
    double norm = 0.;
    for (int d=0; d<rbzone_boundaries_.size2(); ++d) 
      norm += rbzone_boundaries_(i,d)*rbzone_boundaries_(i,d);
    //norm = sqrt(norm);
    double projection = 0.;
    for (int d=0; d<rbzone_boundaries_.size2(); ++d) 
      projection +=  v[d] * rbzone_boundaries_(i,d);
    projection /= norm;
    projection = round(projection);
    if (projection > 1 || projection <= -1) 
      return false;
  }
  return true;
}



bool ClusterHelper::append_column_no_duplicates(Matrix<double>::Type& m, const vector_type& v) const
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



ClusterHelper::vector_type ClusterHelper::map_back_into_brillouin_zone(vector_type v) const
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
  //std::cout << v[0] << " " << v[1] << "\n";
  return v;
}



ClusterHelper::vector_type ClusterHelper::map_back_into_cluster(vector_type v, 
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



unsigned int ClusterHelper::cluster_momentum_approximation(vector_type v) const 
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
    T(0,0)  =  2*M_PI; T(0,1)  =       0; T(0,2)  =       0;
    T(1,0)  =       0; T(1,1)  =  2*M_PI; T(1,2)  =       0;
    T(2,0)  =       0; T(2,1)  =       0; T(2,2)  =  2*M_PI;
    T(3,0)  =  2*M_PI; T(3,1)  =  2*M_PI; T(3,2)  =       0;
    T(4,0)  =  2*M_PI; T(4,1)  =       0; T(4,2)  =  2*M_PI;
    T(5,0)  =       0; T(5,1)  =  2*M_PI; T(5,2)  =  2*M_PI;
    T(6,0)  =  2*M_PI; T(6,1)  = -2*M_PI; T(6,2)  =       0;
    T(7,0)  =  2*M_PI; T(7,1)  =       0; T(7,2)  = -2*M_PI;
    T(8,0)  =       0; T(8,1)  =  2*M_PI; T(8,2)  = -2*M_PI;
    T(9,0)  =  2*M_PI; T(9,1)  =  2*M_PI; T(9,2)  =  2*M_PI;
    T(10,0) = -2*M_PI; T(10,1) =  2*M_PI; T(10,2) =  2*M_PI;
    T(11,0) =  2*M_PI; T(11,1) = -2*M_PI; T(11,2) =  2*M_PI;
    T(12,0) =  2*M_PI; T(12,1) =  2*M_PI; T(12,2) = -2*M_PI;
  }
  for (int i=0; i<cluster_momenta_.size1(); ++i) {
    for (int t=0; t<T.size1(); ++t) {
      for (int s=-1; s<=1; ++s) {
	vector_type diff(dimension());
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



std::vector<unsigned int> ClusterHelper::cluster_neighbors(vector_type k) const 
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
  for (int i=0; i<neighbor.size(); ++i) {
    vector_type nv = matrix_vector_mult(reciprocal_basis_, neighbor[i]);
    nv = map_back_into_brillouin_zone(nv);
    /*if (antiferromagnet) {
      if (!inside_reduced_brillouin_zone(nv)) {
	for (int d=0; d<dimension(); ++d) 
	  nv[d] -= M_PI;
	nv = map_back_into_brillouin_zone(nv);
      }
      if (!inside_reduced_brillouin_zone(nv)) {
	boost::throw_exception(std::logic_error("something doesn't work in ClusterHelper::interpolate)"));
      }
      }*/
    neighbor_index[i] = identify_vertex(nv, cluster_momenta_);
  }
  return neighbor_index;
}



unsigned int ClusterHelper::identify_vertex(const vector_type& v, const Matrix<double>::Type& coordinate_matrix) const
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
  for (int d=0; d<coordinate_matrix.size2(); ++d) 
    sum += coordinate_matrix(i,d);
  sum = round(sum);
  sum *= sum<0 ? -1 : 1;
  return (int(round(sum)))%2;
}


