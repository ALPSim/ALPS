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
  basis_vector_iterator_pair basis_vectors() const {return std::make_pair(basis_.begin(), basis_.end()); }
  basis_vector_iterator_pair reciprocal_basis_vectors() const {return std::make_pair(reciprocal_vectors_.begin(), reciprocal_vectors_.end()); }
  basis_vector_iterator_pair cell_origins() const{return std::make_pair(cell_origins_.begin(), cell_origins_.end());}
  offset_iterator_pair cell_offsets() const{return std::make_pair(cell_offsets_.begin(), cell_offsets_.end());}
  momenta_iterator_pair cell_momenta() const{return std::make_pair(momenta_.begin(), momenta_.end());}
  void compute_unit_cell_volume(){
    if(dimension_==1){ unit_cell_volume_=basis_[0][0]; }
    else if(dimension_==2){
      unit_cell_volume_=std::abs(basis_[0][0]*basis_[1][1]-basis_[0][1]*basis_[1][0]);
    }else if(dimension_==3){
      unit_cell_volume_=std::abs(basis_[2][0]*(basis_[0][1]*basis_[1][2]-basis_[1][1]*basis_[0][2])+basis_[2][1]*(basis_[0][2]*basis_[1][0]-basis_[0][0]*basis_[1][2])+basis_[2][2]*(basis_[0][0]*basis_[1][1]-basis_[0][1]*basis_[1][0]));
    }
  }
private:
  int volume_;
  int unit_cell_volume_;
  int dimension_;
  std::vector<std::vector<double> > basis_;
  std::vector<std::vector<double> > reciprocal_vectors_;
  std::vector<std::vector<double> > cell_origins_;
  offset_vector_type cell_offsets_;
  momenta_vector_type momenta_;
};
