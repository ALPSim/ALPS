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

  alps_graph(const alps::params& params):graph_(alps::make_deprecated_parameters(params)){
    //set up cell origins and offsets
    cell_iterator itc, itc_end;
    for (boost::tie(itc, itc_end) = graph_.cells(); itc != itc_end; ++itc) {
      cell_origins_.push_back(graph_.origin(*itc));
      cell_offsets_.push_back(graph_.offset(*itc));
    }
    //set up momenta
    momentum_iterator itm, itm_end;
    for (boost::tie(itm, itm_end) = graph_.momenta(); itm != itm_end; ++itm) {
      momenta_.push_back(graph_.momentum(*itm));
    }
    //buffer volume and dimension
    volume_=graph_.volume();
    dimension_=graph_.dimension();
    
    //buffer the basis vectors
    basis_vector_iterator_pair BVIP=graph_.basis_vectors();
    for(basis_vector_iterator it=BVIP.first; it!=BVIP.second;++it){
      basis_.push_back(*it);
    }
    basis_vector_iterator_pair RVIP=graph_.reciprocal_basis_vectors();
    for(basis_vector_iterator it=RVIP.first; it!=RVIP.second;++it){
      reciprocal_vectors_.push_back(*it);
    }

  }
  virtual ~alps_graph(){}
   int volume() const{return volume_;}
   int dimension() const{return dimension_;}
   basis_vector_iterator_pair basis_vectors() const {return std::make_pair(basis_.begin(), basis_.end()); }
   basis_vector_iterator_pair reciprocal_basis_vectors() const {return std::make_pair(reciprocal_vectors_.begin(), reciprocal_vectors_.end()); }
   basis_vector_iterator_pair cell_origins() const{return std::make_pair(cell_origins_.begin(), cell_origins_.end());}
   offset_iterator_pair cell_offsets() const{return std::make_pair(cell_offsets_.begin(), cell_offsets_.end());}
   momenta_iterator_pair cell_momenta() const{return std::make_pair(momenta_.begin(), momenta_.end());}
private:
  int volume_;
  int dimension_;
  const alps::graph_helper<alps::graph_helper<>::graph_type> graph_;
  std::vector<std::vector<double> > basis_;
  std::vector<std::vector<double> > reciprocal_vectors_;
  std::vector<std::vector<double> > cell_origins_;
  offset_vector_type cell_offsets_;
  momenta_vector_type momenta_;
};
