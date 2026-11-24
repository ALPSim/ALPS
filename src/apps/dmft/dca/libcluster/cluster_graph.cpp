/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2026 by Emanuel Gull
 *
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/


#include<sstream>
#include<cctype>
#include<cmath>
#include<climits>
#include "cluster_graph.h"
#include <iostream>

namespace {
// The value of the attribute `name` in the start tag `tag`, e.g. name="..." or
// name = "...", or "" if the tag has no such attribute.
std::string attribute(const std::string &tag, const std::string &name) {
  std::size_t pos = 0;
  while ((pos = tag.find(name, pos)) != std::string::npos) {
    const bool starts_word = pos > 0 && std::isspace(static_cast<unsigned char>(tag[pos - 1]));
    std::size_t p = pos + name.size();
    while (p < tag.size() && std::isspace(static_cast<unsigned char>(tag[p]))) ++p;
    if (starts_word && p < tag.size() && tag[p] == '=') {
      ++p;
      while (p < tag.size() && std::isspace(static_cast<unsigned char>(tag[p]))) ++p;
      if (p < tag.size() && (tag[p] == '"' || tag[p] == '\'')) {
        const std::size_t end = tag.find(tag[p], p + 1);
        if (end != std::string::npos) return tag.substr(p + 1, end - p - 1);
      }
    }
    pos += name.size();
  }
  return "";
}
}

// Reads the basis of the superlattice "<cluster_name> superlattice" from a
// cluster library in the ALPS lattice library format:
//   <LATTICE name="Betts2D-8A superlattice" dimension="2">
//     <BASIS> <VECTOR>2 2</VECTOR> <VECTOR>2 -2</VECTOR> </BASIS> ...
// The name has to match exactly, and each basis vector needs `dimension`
// integer entries (#443).
std::vector<std::vector<double> > alps_graph::lookup_basis_from_name(const std::string &cluster_name, const std::string &cluster_library) const{
  std::ifstream xml_file(cluster_library.c_str());
  if(!xml_file.good()) throw std::runtime_error("cluster library: "+cluster_library+" could not be opened.");
  std::stringstream contents;
  contents << xml_file.rdbuf();
  std::string xml = contents.str();
  // Entries inside <!-- --> are not read: blank out the comments, keeping
  // the positions of everything else. An unterminated comment runs to the end.
  for (std::size_t begin = xml.find("<!--"); begin != std::string::npos; begin = xml.find("<!--", begin)) {
    const std::size_t close = xml.find("-->", begin + 4);
    const std::size_t end = close == std::string::npos ? xml.size() : close + 3;
    xml.replace(begin, end - begin, end - begin, ' ');
  }
  const std::string where = "cluster " + cluster_name + " in " + cluster_library + ": ";
  const std::string wanted = cluster_name + " superlattice";

  std::size_t found = std::string::npos, found_end = std::string::npos;
  for (std::size_t pos = xml.find("<LATTICE"); pos != std::string::npos; pos = xml.find("<LATTICE", pos + 1)) {
    const std::size_t after = pos + 8;
    if (after >= xml.size() || !(std::isspace(static_cast<unsigned char>(xml[after])) || xml[after] == '>'))
      continue;  // <LATTICEGRAPH ...>
    const std::size_t tag_end = xml.find('>', pos);
    if (tag_end == std::string::npos) break;
    if (attribute(xml.substr(pos, tag_end - pos), "name") != wanted) continue;
    if (found != std::string::npos)
      throw std::runtime_error(where + "the library defines it more than once.");
    found = pos;
    found_end = tag_end;
  }
  if (found != std::string::npos && xml[found_end - 1] == '/') {
    // <LATTICE .../> has no basis; the next lattice's must not be taken for it
    throw std::runtime_error(where + "has no <BASIS> ... </BASIS>.");
  }
  if (found == std::string::npos)
    throw std::runtime_error(where + "not found. dca.LATTICE has to be the name of a LATTICE \""
                             + wanted + "\" in the library, exactly.");

  const std::string dim_string = attribute(xml.substr(found, found_end - found), "dimension");
  if (dim_string != "1" && dim_string != "2" && dim_string != "3")
    throw std::runtime_error(where + "dimension=\"" + dim_string + "\"; clusters have dimension 1, 2 or 3.");
  const int dimension = dim_string[0] - '0';

  const std::size_t lattice_end = xml.find("</LATTICE>", found_end);
  const std::size_t basis_begin = xml.find("<BASIS>", found_end);
  const std::size_t basis_end = xml.find("</BASIS>", found_end);
  if (lattice_end == std::string::npos || basis_begin == std::string::npos || basis_end == std::string::npos
      || basis_begin > lattice_end || basis_end > lattice_end || basis_end < basis_begin)
    throw std::runtime_error(where + "has no <BASIS> ... </BASIS>.");

  std::vector<std::vector<double> > basis;
  for (std::size_t pos = xml.find("<VECTOR>", basis_begin); pos != std::string::npos && pos < basis_end;
       pos = xml.find("<VECTOR>", pos + 1)) {
    const std::size_t end = xml.find("</VECTOR>", pos);
    if (end == std::string::npos || end > basis_end)
      throw std::runtime_error(where + "a <VECTOR> of the basis has no </VECTOR>.");
    const std::string text = xml.substr(pos + 8, end - pos - 8);
    std::stringstream entries(text);
    std::vector<double> vector;
    std::string entry;
    while (entries >> entry) {
      std::size_t used = 0;
      double value = 0;
      try { value = std::stod(entry, &used); } catch (std::exception const &) { used = 0; }
      if (used != entry.size() || !std::isfinite(value) || value != std::round(value) || std::abs(value) > INT_MAX)
        throw std::runtime_error(where + "basis vector \"" + text + "\": \"" + entry + "\" is not an integer."
                                 + " dca needs a numeric basis of integer lattice vectors"
                                 + " (the 1D chains are Chain1D_<N> in cluster1d.xml).");
      vector.push_back(value);
    }
    if (int(vector.size()) != dimension)
      throw std::runtime_error(where + "basis vector \"" + text + "\" has " + std::to_string(vector.size())
                               + " entries, dimension=\"" + dim_string + "\" needs " + dim_string + ".");
    basis.push_back(vector);
  }
  if (int(basis.size()) != dimension)
    throw std::runtime_error(where + "the basis has " + std::to_string(basis.size()) + " vectors, dimension=\""
                             + dim_string + "\" needs " + dim_string + ".");
  return basis;
}
void alps_graph::compute_unit_cell_volume(){
  if(dimension_==1){ unit_cell_volume_= std::abs(basis_vectors_[0][0]); }
  else if(dimension_==2){
    unit_cell_volume_=std::abs(basis_vectors_[0][0]* basis_vectors_[1][1]- basis_vectors_[0][1]* basis_vectors_[1][0]);
  }else if(dimension_==3){
    std::vector<double> c=cross(basis_vectors_[0], basis_vectors_[1]);
    unit_cell_volume_=std::abs(basis_vectors_[2][0]*c[0]+ basis_vectors_[2][1]*c[1]+ basis_vectors_[2][2]*c[2]);
  }
}
std::vector<std::vector<double> > alps_graph::compute_reciprocal_from_basis(const std::vector<std::vector<double> > &basis) const{
  std::vector<std::vector<double> > basis_in_3d(3, std::vector<double>(3, 0.));
  std::vector<std::vector<double> > reciprocal_basis_in_3d(3, std::vector<double>(3, 0.));
  std::vector<std::vector<double> > reciprocal_basis(dimension(), std::vector<double>(dimension(), 0.));
  
  for(int i=0;i<dimension();++i){
    for(int j=0;j<dimension();++j){
      basis_in_3d[i][j]=basis[i][j];
    }
  }for(int i=dimension();i<3;++i){
    basis_in_3d[i][i]=1.;
  }
  std::vector<double> c0=cross(basis_in_3d[1], basis_in_3d[2]);
  std::vector<double> c1=cross(basis_in_3d[2], basis_in_3d[0]);
  std::vector<double> c2=cross(basis_in_3d[0], basis_in_3d[1]);
  double numerator=c0[0]*basis_in_3d[0][0]+c0[1]*basis_in_3d[0][1]+c0[2]*basis_in_3d[0][2];

  reciprocal_basis_in_3d[0][0]=2*M_PI/numerator*c0[0];
  reciprocal_basis_in_3d[0][1]=2*M_PI/numerator*c0[1];
  reciprocal_basis_in_3d[0][2]=2*M_PI/numerator*c0[2];
  reciprocal_basis_in_3d[1][0]=2*M_PI/numerator*c1[0];
  reciprocal_basis_in_3d[1][1]=2*M_PI/numerator*c1[1];
  reciprocal_basis_in_3d[1][2]=2*M_PI/numerator*c1[2];
  reciprocal_basis_in_3d[2][0]=2*M_PI/numerator*c2[0];
  reciprocal_basis_in_3d[2][1]=2*M_PI/numerator*c2[1];
  reciprocal_basis_in_3d[2][2]=2*M_PI/numerator*c2[2];
  
  for(int i=0;i<dimension();++i){
    for(int j=0;j<dimension();++j){
      reciprocal_basis[i][j]=reciprocal_basis_in_3d[i][j];
    }
  }
  return reciprocal_basis;
}

namespace {
// Equal sizes, and entries within 1e-8.
template <typename T>
bool same_table(const std::vector<std::vector<T> > &a, const std::vector<std::vector<T> > &b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (a[i].size() != b[i].size()) return false;
    for (std::size_t j = 0; j < a[i].size(); ++j)
      if (std::abs(double(a[i][j]) - double(b[i][j])) > 1.e-8) return false;
  }
  return true;
}
}

bool alps_graph::operator==(const alps_graph &g) const{
  return volume_ == g.volume_ && unit_cell_volume_ == g.unit_cell_volume_ && dimension_ == g.dimension_
      && same_table(basis_vectors_, g.basis_vectors_)
      && same_table(reciprocal_basis_vectors_, g.reciprocal_basis_vectors_)
      && same_table(cell_origins_, g.cell_origins_)
      && same_table(cell_offsets_, g.cell_offsets_)
      && same_table(momenta_, g.momenta_);
}

void alps_graph::build(const std::vector<std::vector<double> > basis, int L){
  basis_vectors_ =basis;
  dimension_=basis.size(); if(dimension_>3 || dimension_<1) throw std::invalid_argument("clusters work for dimensions 1, 2, and 3.");
  volume_=dimension_==1?L:dimension_==2?L*L:L*L*L;
  compute_unit_cell_volume();
  if(unit_cell_volume_==0)
    throw std::invalid_argument("cluster: the basis vectors are linearly dependent, the cluster has no sites.");
  
  reciprocal_basis_vectors_ =compute_reciprocal_from_basis(basis);
  std::pair<std::vector<std::vector<double> >,std::vector<std::vector<int> > > origins_and_offsets=compute_cell_origins_and_offsets(L);
  cell_origins_=origins_and_offsets.first;
  cell_offsets_=origins_and_offsets.second;
  momenta_=compute_momenta(L);
}

std::vector<std::vector<double> > alps_graph::compute_momenta(int L)const{
  std::vector<std::vector<double> > momenta;
  std::vector<double> momentum(dimension(), 0.);
  if(dimension()==3){
    for(int i=0;i<L;++i){
      for(int j=0;j<L;++j){
        for(int k=0;k<L;++k){
          int offset_x=i>L/2?i-L:i;
          int offset_y=j>L/2?j-L:j;
          int offset_z=k>L/2?k-L:k;
          momentum[0]=
                  reciprocal_basis_vectors_[0][0]*offset_x/(double)L+ reciprocal_basis_vectors_[1][0]*offset_y/(double)L+
                  reciprocal_basis_vectors_[2][0]*offset_z/(double)L;
          momentum[1]=
                  reciprocal_basis_vectors_[0][1]*offset_x/(double)L+ reciprocal_basis_vectors_[1][1]*offset_y/(double)L+
                  reciprocal_basis_vectors_[2][1]*offset_z/(double)L;
          momentum[2]=
                  reciprocal_basis_vectors_[0][2]*offset_x/(double)L+ reciprocal_basis_vectors_[1][2]*offset_y/(double)L+
                  reciprocal_basis_vectors_[2][2]*offset_z/(double)L;
          momenta.push_back(momentum);
        }
      }
    }
  }else if(dimension()==2){
    for(int i=0;i<L;++i){
      for(int j=0;j<L;++j){
        int offset_x=i>L/2?i-L:i;
        int offset_y=j>L/2?j-L:j;
        momentum[0]=
                reciprocal_basis_vectors_[0][0]*offset_x/(double)L+ reciprocal_basis_vectors_[1][0]*offset_y/(double)L;
        momentum[1]=
                reciprocal_basis_vectors_[0][1]*offset_x/(double)L+ reciprocal_basis_vectors_[1][1]*offset_y/(double)L;
        momenta.push_back(momentum);
      }
    }
  }else if(dimension()==1){
    for(int i=0;i<L;++i){
      int offset_x=i>L/2?i-L:i;
      momentum[0]= reciprocal_basis_vectors_[0][0]*offset_x/(double)L;
      momenta.push_back(momentum);
    }
  }
  return momenta;
}
std::pair<std::vector<std::vector<double> >,std::vector<std::vector<int> > > alps_graph::compute_cell_origins_and_offsets(int L){
  int ncells=dimension()==1?L:dimension()==2?L*L:dimension()==3?L*L*L:throw std::runtime_error("dimension() should be 1, 2, or 3");
  std::vector<std::vector<double> > cell_origins(ncells, std::vector<double>(dimension(), 0.));
  std::vector<std::vector<int   > > cell_offsets(ncells, std::vector<int>(dimension(), 0));
  std::vector<double> cell_origin(dimension(), 0.);
  std::vector<int   > cell_offset(dimension(), 0 );
  if(dimension()==3){
    for(int i=0;i<L;++i){
      for(int j=0;j<L;++j){
        for(int k=0;k<L;++k){
          cell_origin[0] =i* basis_vectors_[0][0]+j* basis_vectors_[1][0]+k* basis_vectors_[2][0];
          cell_origin[1] =i* basis_vectors_[0][1]+j* basis_vectors_[1][1]+k* basis_vectors_[2][1];
          cell_origin[2] =i* basis_vectors_[0][2]+j* basis_vectors_[1][2]+k* basis_vectors_[2][2];
          cell_offset[0] =i; cell_offset[1] =j; cell_offset[2] =k;
          cell_origins[i*L*L+j*L+k]=cell_origin;
          cell_offsets[i*L*L+j*L+k]=cell_offset;
        }
      }
    }
  }else if(dimension()==2){
    for(int i=0;i<L;++i){
      for(int j=0;j<L;++j){
        cell_origin[0] =i* basis_vectors_[0][0]+j* basis_vectors_[1][0];
        cell_origin[1] =i* basis_vectors_[0][1]+j* basis_vectors_[1][1];
        cell_offset[0] =i; cell_offset[1] =j;
        cell_origins[i*L+j]=cell_origin;
        cell_offsets[i*L+j]=cell_offset;
      }
    }
  }else if(dimension()==1){
    for(int i=0;i<L;++i){
      cell_origin[0] =i* basis_vectors_[0][0];
      cell_offset[0] =i;
      cell_origins[i]=cell_origin;
      cell_offsets[i]=cell_offset;
    }
  }
  return std::make_pair(cell_origins, cell_offsets);
}

