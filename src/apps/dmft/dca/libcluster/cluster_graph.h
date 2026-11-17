/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2015 by Emanuel Gull <egull@umich.edu>
 *
 *
 *****************************************************************************/

#pragma once
#include<vector>
#include<string>
#include<fstream>
#include<stdexcept>

class alps_graph {
public:
  //to disentangle things from ALPS we have to hide the cell iterator, offset type, and momentum iterator calls.
  typedef std::vector<std::vector<double> > basis_vector;
  typedef basis_vector::const_iterator basis_vector_iterator;
  typedef std::pair<basis_vector_iterator,basis_vector_iterator> basis_vector_iterator_pair;
  
  typedef std::vector<std::vector<int> > offset_vector;
  typedef offset_vector::const_iterator offset_vector_iterator;
  typedef std::pair<offset_vector_iterator,offset_vector_iterator> offset_vector_iterator_pair;

  ///Constructor using a name for the cluster, an xml path, and a lattice extent L.
  alps_graph(const std::string &cluster_name, const std::string &cluster_library, int L){
    build(lookup_basis_from_name(cluster_name, cluster_library), L);
  }
  ///return the volume of the lattice. this is L, L*L, or L*L*L, depending on the dimension.
  int volume() const{return volume_;}
  ///dimension. Should be 1, 2, or 3.
  int dimension() const{return dimension_;}
  ///volume of the unit cell (aka the 'cluster size', i.e. 8 for an 8-site cluster)
  int unit_cell_volume() const { return unit_cell_volume_;}
  ///const reference to the basis vectors
  const basis_vector &basis_vectors() const {return basis_vectors_; }
  ///const reference to the basis vectors in reciprocal space
  const basis_vector &reciprocal_basis_vectors() const {return reciprocal_basis_vectors_; }
  ///returns iterators to the beginning and end of the cell origin vectors.
  const basis_vector &cell_origins() const{return cell_origins_;}
  ///returns iterators to the beginning and end of the cell offset vectors.
  const offset_vector &cell_offsets() const{return cell_offsets_;}
  ///returns iterators to the beginning and end of the cell momenta vectors. TODO: remove this.
  const basis_vector &cell_momenta() const{return momenta_;}
  ///compute and set the unit cell volume ('8' for an 8-site cluster)
  void compute_unit_cell_volume();
  ///given a basis, compute its reciprocal basis. This works for d=1, 2, and 3.
  std::vector<std::vector<double> > compute_reciprocal_from_basis(const std::vector<std::vector<double> > &basis) const;
  ///compute the cell origins (i.e. the (x,y,z) positions of the translated unit cells and their translation offsets
  std::pair<std::vector<std::vector<double> >,std::vector<std::vector<int> > > compute_cell_origins_and_offsets(int L);
  ///compute all possible momenta explicitly (rahter than via ALPS)
  std::vector<std::vector<double> > compute_momenta(int L)const;
  ///comparison operator, mainly used for tests.
  bool operator==(const alps_graph &g) const;
  ///if the reciprocal basis is correct, this test will validate that each vector is orthogonal to all other basis vectors
  double validate_orthogonality() const;
private:
  ///the cross product in 3d.
  std::vector<double> cross(const std::vector<double> &v1, const std::vector<double> &v2) const{
    std::vector<double> c(3, 0.);
    c[0]=v1[1]*v2[2]-v1[2]*v2[1];
    c[1]=v1[2]*v2[0]-v1[0]*v2[2];
    c[2]=v1[0]*v2[1]-v1[1]*v2[0];
    return c;
  }
  ///given a basis and a number of unit cells L, build a lattice.
  void build(const std::vector<std::vector<double> > basis, int L);
  ///this function builds a cluster based on an xml cluster description file in the ALPS lattice library format.
  std::vector<std::vector<double> > lookup_basis_from_name(const std::string &cluster_name, const std::string &cluster_library) const;
  ///volume of the lattice (# of unit cells)
  int volume_;
  ///volume of a single unit cell
  int unit_cell_volume_;
  ///dimension, is 1, 2, or 3.
  int dimension_;
  ///a vector of basis vectors
  std::vector<std::vector<double> > basis_vectors_;
  ///a vector of reciprocal basic vectors in the Brillouin zone
  std::vector<std::vector<double> > reciprocal_basis_vectors_;
  ///a vector of locations of the origin of translated unit cells
  std::vector<std::vector<double> > cell_origins_;
  ///a vector of translation vectors of the translated unit cells
  std::vector<std::vector<int   > > cell_offsets_;
  ///points in the brillouin zone of the cluster, which is smaller than the BZ of the lattice (-pi to pi) but tiles it when shifted by reciprocal basis vectors
  std::vector<std::vector<double> > momenta_;
};
