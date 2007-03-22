/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2007 by Sebastian Fuchs <sfuchs@gwdg.de>
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

#ifndef ALPS_CLUSTERTRANSFORMER_H
#define ALPS_CLUSTERTRANSFORMER_H

#include "types.h"
#include "green_function.h"
#include "hilberttransformer.h"
#include <alps/lattice/graph_helper.h>



class ClusterHelper: public alps::graph_helper<alps::graph_helper<>::graph_type>
{
public:

  ClusterHelper(const alps::Parameters& params);

  matsubara_green_function_t transform_into_real_space(const matsubara_green_function_t& G) const;
  itime_green_function_t transform_into_real_space(const itime_green_function_t& G) const;
  matsubara_green_function_t transform_into_k_space(const matsubara_green_function_t& G) const;
  itime_green_function_t transform_into_k_space(const itime_green_function_t& G) const;
  
protected:
  
  std::pair<vertex_iterator, vertex_iterator> cluster_sites() const{
    return alps::detail::vertices_wrap(unit_cell().graph());
  }
  
  vector_type cluster_coordinate(const vertex_descriptor itv) const{
    return boost::get(alps::coordinate_t(), unit_cell().graph(), itv);
  }

  vector_type cluster_momentum(const vertex_descriptor itv) const{
    vector_type v = cluster_coordinate(itv);
    vector_type k(dimension());
    basis_vector_iterator itb, itb_end;
    boost::tie(itb, itb_end) = basis_vectors();
    for (uint i=0; itb!=itb_end; ++itb, ++i) {
      double norm = 0.;
      for(uint j=0; j<dimension(); ++j) 
	norm += (*itb)[j] * (*itb)[j];
      k[i] = 0.;
      for(uint j=0; j<dimension(); ++j) {
	k[i] += 2*M_PI*(*itb)[j]*v[j]/norm;
      }
    }
    return k;
  }

uint n_cluster_sites() const { 
  return (num_sites() / volume()); }

};


class ClusterTransformer : public ClusterHelper
{
public:
  ClusterTransformer(const alps::Parameters& params) :
    ClusterHelper(params) {}

  virtual matsubara_green_function_t operator()(matsubara_green_function_t & G_omega, 
						matsubara_green_function_t &G0_omega, 
						const double mu, const double h, 
						const double beta,const bool paramagnet)=0; 

  matsubara_green_function_t calculate_initial_G0(const alps::Parameters& parms);
};


class DCATransformer: public ClusterTransformer
{
public:

  DCATransformer(const alps::Parameters& params) :
    ClusterTransformer(params), t_(boost::lexical_cast<double>(params["t"])) {}


  matsubara_green_function_t operator()(matsubara_green_function_t &G_omega, 
					matsubara_green_function_t &G0_omega, 
					const double mu, const double h, const double beta, const bool paramagnet);

private:
  double epsilon(const vector_type k) const
  {
    double eps = 0;
    for (uint i=0; i<dimension(); ++i) 
      eps -= 2*t_*cos (k[i]);
    return eps; 
  }
  
  const double t_;
};


#endif /*ALPS_CLUSTERTRANSFORMER_H*/
