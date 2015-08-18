/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2015 by Emanuel Gull <egull@umich.edu>
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

#pragma once
#include "cluster.h"
#include "alps/lattice/graph_helper.h"
#include "alps/numeric/vector_functions.hpp"
#include <alps/ngs/make_deprecated_parameters.hpp>

class alps_graph:public abstract_graph{
public:
  //to disentangle things from ALPS we have to hide the cell iterator, offset type, and momentum iterator calls.
  typedef alps::graph_helper<alps::graph_helper<>::graph_type>::cell_iterator cell_iterator;
  typedef std::pair<cell_iterator,cell_iterator> cell_iterator_pair;
  typedef alps::graph_helper<alps::graph_helper<>::graph_type>::cell_descriptor cell_descriptor;
  typedef alps::graph_helper<alps::graph_helper<>::graph_type>::momentum_iterator momentum_iterator;
  typedef std::pair<momentum_iterator,momentum_iterator> momentum_iterator_pair;
  typedef std::vector<std::vector<double> > basis_vector;
  typedef basis_vector::const_iterator basis_vector_iterator;
  typedef std::pair<basis_vector_iterator,basis_vector_iterator> basis_vector_iterator_pair;
  
  alps_graph(const alps::params& params){
    const alps::graph_helper<alps::graph_helper<>::graph_type> graph(alps::make_deprecated_parameters(params));
    //set up cell origins and offsets
    cell_iterator itc, itc_end;
    for (boost::tie(itc, itc_end) = graph.cells(); itc != itc_end; ++itc) {
      cell_origins_.push_back(graph.origin(*itc));
      cell_offsets_.push_back(graph.offset(*itc));
    }
    //set up momenta
    momentum_iterator itm, itm_end;
    for (boost::tie(itm, itm_end) = graph.momenta(); itm != itm_end; ++itm) {
      momenta_.push_back(graph.momentum(*itm));
    }
    //buffer volume and dimension
    volume_=graph.volume();
    dimension_=graph.dimension();
    if(graph.dimension()>3) throw std::invalid_argument("we only deal with graphs of dimension 1, 2, or 3");
    
    //buffer the basis vectors
    basis_vector_iterator_pair BVIP=graph.basis_vectors();
    for(basis_vector_iterator it=BVIP.first; it!=BVIP.second;++it){
      basis_.push_back(*it);
    }
    basis_vector_iterator_pair RVIP=graph.reciprocal_basis_vectors();
    for(basis_vector_iterator it=RVIP.first; it!=RVIP.second;++it){
      reciprocal_vectors_.push_back(*it);
    }
    //compute the unit cell volume
    compute_unit_cell_volume();
  }
  virtual ~alps_graph(){}
  int volume() const{return volume_;}
  int dimension() const{return dimension_;}
  int unit_cell_volume(){ return unit_cell_volume_;}
  //if the reciprocal basis is correct, this test will validate that each vector is orthogonal to all other basis vectors
  double validate_orghogonality(){
    double deviation_from_orthogonality=0.;
    for(int i=0;i<basis_.size();++i){
      for(int j=0;j<reciprocal_vectors_.size();++j){
        double inner_product=0.;
        for(int k=0;k<dimension();++k){
          inner_product+=basis_[i][k]*reciprocal_vectors_[j][k];
        }
        deviation_from_orthogonality+=inner_product-2.*M_PI*(i==j);
      }
    }
    return deviation_from_orthogonality;
  }
  basis_vector_iterator_pair basis_vector_iterators() const {return std::make_pair(basis_.begin(), basis_.end()); }
  const basis_vector &basis_vectors() const {return basis_; }
  const basis_vector &reciprocal_basis_vectors() const {return reciprocal_vectors_; }
  basis_vector_iterator_pair reciprocal_basis_vector_iterators() const {return std::make_pair(reciprocal_vectors_.begin(), reciprocal_vectors_.end()); }
  basis_vector_iterator_pair cell_origins() const{return std::make_pair(cell_origins_.begin(), cell_origins_.end());}
  offset_iterator_pair cell_offsets() const{return std::make_pair(cell_offsets_.begin(), cell_offsets_.end());}
  momenta_iterator_pair cell_momenta() const{return std::make_pair(momenta_.begin(), momenta_.end());}
  void compute_unit_cell_volume(){
    if(dimension_==1){ unit_cell_volume_=basis_[0][0]; }
    else if(dimension_==2){
      unit_cell_volume_=std::abs(basis_[0][0]*basis_[1][1]-basis_[0][1]*basis_[1][0]);
    }else if(dimension_==3){
      std::vector<double> c=cross(basis_[0], basis_[1]);
      unit_cell_volume_=std::abs(basis_[2][0]*c[0]+basis_[2][1]*c[1]+basis_[2][2]*c[2]);
    }
  }
  std::vector<std::vector<double> > compute_reciprocal_from_basis(const std::vector<std::vector<double> > &basis){
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
    using alps::numeric::operator*;
    reciprocal_basis_in_3d[0]=2*M_PI/numerator*c0;
    reciprocal_basis_in_3d[1]=2*M_PI/numerator*c1;
    reciprocal_basis_in_3d[2]=2*M_PI/numerator*c2;

    for(int i=0;i<dimension();++i){
      for(int j=0;j<dimension();++j){
        reciprocal_basis[i][j]=reciprocal_basis_in_3d[i][j];
      }
    }
    return reciprocal_basis;
  }
  std::pair<std::vector<std::vector<double> >,std::vector<std::vector<double> > > cell_origins_and_offsets_by_hand(int L){
    std::vector<std::vector<double> > cell_origins;
    std::vector<std::vector<double> > cell_offsets;
    std::vector<double> cell_origin(dimension(), 0.);
    std::vector<double> cell_offset(dimension(), 0.);
    if(dimension()==3){
      for(int i=0;i<L;++i){
        for(int j=0;j<L;++j){
          for(int k=0;k<L;++k){
            cell_origin[0] =i*basis_[0][0]+j*basis_[1][0]+k*basis_[2][0];
            cell_origin[1] =i*basis_[0][1]+j*basis_[1][1]+k*basis_[2][1];
            cell_origin[2] =i*basis_[0][2]+j*basis_[1][2]+k*basis_[2][2];
            cell_offset[0] =i; cell_offset[1] =j; cell_offset[2] =k;
            cell_origins.push_back(cell_origin);
            cell_offsets.push_back(cell_offset);
          }
        }
      }
    }else if(dimension()==2){
      for(int i=0;i<L;++i){
        for(int j=0;j<L;++j){
          cell_origin[0] =i*basis_[0][0]+j*basis_[1][0];
          cell_origin[1] =i*basis_[0][1]+j*basis_[1][1];
          cell_offset[0] =i; cell_offset[1] =j;
          cell_origins.push_back(cell_origin);
          cell_offsets.push_back(cell_offset);
        }
      }
    }else if(dimension()==1){
      for(int i=0;i<L;++i){
        cell_origin[0] =i*basis_[0][0];
        cell_offset[0] =i;
        cell_origins.push_back(cell_origin);
        cell_offsets.push_back(cell_offset);
      }
    }
    return std::make_pair(cell_origins, cell_offsets);
  }
private:
  std::vector<double> cross(const std::vector<double> &v1, const std::vector<double> &v2){
    std::vector<double> c(3, 0.);
    c[0]=v1[1]*v2[2]-v1[2]*v2[1];
    c[1]=v1[2]*v2[0]-v1[0]*v2[2];
    c[2]=v1[0]*v2[1]-v1[1]*v2[0];
    return c;
  }
  int volume_;
  int unit_cell_volume_;
  int dimension_;
  std::vector<std::vector<double> > basis_;
  std::vector<std::vector<double> > reciprocal_vectors_;
  std::vector<std::vector<double> > cell_origins_;
  std::vector<std::vector<int> >    cell_offsets_;
  momenta_vector_type momenta_;
};
