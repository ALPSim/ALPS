/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2006 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabian Stoeckli <fabstoec@student.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_APPLICATIONS_MC_SPIN_LOCALUPDATE_H_
#define ALPS_APPLICATIONS_MC_SPIN_LOCALUPDATE_H_

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/random/uniform_int.hpp>
#include <vector>

#include "potts.h"
#include "ising.h"
#include "on.h"
#include "matrices.h"

// Diagonal and single-element couplings are symmetric.
template <class MAT>
const MAT& transposed(const MAT& J) { return J; }

template <class T, int N>
SquareMatrix<T,N> transposed(const SquareMatrix<T,N>& J) {
  SquareMatrix<T,N> result;
  for (int i=0; i<N; ++i)
    for (int j=0; j<N; ++j)
      result.setElement_nc(i,j,J[j*N+i]);
  return result;
}

// The coupling of an out-edge, as a matrix whose row index belongs to vertex.
// J couples the bond's stored source to its stored target, but an out-edge of
// an undirected graph always has vertex as its source. bond_source holds the
// stored source of each bond, indexed by edge_index.
template <class Graph, class CouplingMap>
typename CouplingMap::value_type
oriented_coupling(const Graph& graph, const CouplingMap& coupling,
                  const std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>& bond_source,
                  typename boost::graph_traits<Graph>::vertex_descriptor vertex,
                  typename boost::graph_traits<Graph>::edge_descriptor edge)
{
  typename CouplingMap::value_type J = coupling[edge];
  if (bond_source[boost::get(boost::edge_index,graph,edge)] != vertex)
    J = transposed(J);
  return J;
}

template <class Graph, class MomentMap, class RNG, class CouplingMap, class SpinfactorMap>
void local_update(const Graph& graph, MomentMap& moment,
                  double beta, RNG& rng,  const CouplingMap& coupling,
                  const std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>& bond_source,
                  const SpinfactorMap& spinfactor, double g,
                  TinyVector<double,CouplingMap::value_type::dim>
                    h = TinyVector<double,CouplingMap::value_type::dim>(0.0))

{
  // some typedefs for clarity
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_traits<MomentMap>::value_type moment_type;
  typedef typename moment_type::update_type update_type;
  typedef typename CouplingMap::value_type MAT;

  // a local update consists of on averaqe one update per site
  for (int steps=0; steps<boost::num_vertices(graph);++steps) {
    // step i)
    vertex_descriptor vertex=*(boost::vertices(graph).first 
                              + int (rng()*boost::num_vertices(graph)));

    // step ii)
    update_type update=moment[vertex].random_update(rng);

    // step iii)
    double delta_E=0;
    typename boost::graph_traits<Graph>::out_edge_iterator edge,end;
    for (boost::tie(edge,end)=boost::out_edges(vertex,graph);edge!=end;++edge) {
      double fact = -spinfactor[boost::source(*edge,graph)]*spinfactor[boost::target(*edge,graph)];
      MAT J = oriented_coupling(graph,coupling,bond_source,vertex,*edge);
      J = J * fact;
          delta_E += bond_energy_change(moment[boost::source(*edge,graph)],
                               moment[boost::target(*edge,graph)],J,update);
    }
    delta_E += spinfactor[vertex]*site_energy_change(moment[vertex], h, update)  * g;
    
    // step iv)
    if (delta_E<0 || rng() < exp(-beta*delta_E)) 
      moment[vertex].update(update);
  }
}

// generalized version, with onsite_energy_terms
template <class Graph, class MomentMap, class RNG, class CouplingMap, 
            class SelfinteractionMap, class SpinfactorMap>
void local_update_self(const Graph& graph, MomentMap& moment,
                  double beta, RNG& rng,  const CouplingMap& coupling, 
                  const std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>& bond_source,
                  const SelfinteractionMap& selfinteraction, 
                  const SpinfactorMap& spinfactor, double g,
                  TinyVector<double,CouplingMap::value_type::dim> 
                    h = TinyVector<double,CouplingMap::value_type::dim>(0.0))
{
  // some typedefs for clarity
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_traits<MomentMap>::value_type moment_type;
  typedef typename moment_type::update_type update_type;
  typedef typename CouplingMap::value_type MAT;

  // a local update consists of on averaqe one update per site
  for (int steps=0; steps<boost::num_vertices(graph);++steps) {
    // step i)
    vertex_descriptor vertex=*(boost::vertices(graph).first + 
                                 int(rng()*boost::num_vertices(graph)));
    
    // step ii)
    update_type update = moment[vertex].random_update(rng);
  
    // step iii)
    double delta_E=0;
    typename boost::graph_traits<Graph>::out_edge_iterator edge,end;
    for (boost::tie(edge,end)=boost::out_edges(vertex,graph);edge!=end;++edge) {
      double fact = -spinfactor[boost::source(*edge,graph)]*spinfactor[boost::target(*edge,graph)];
      MAT J = oriented_coupling(graph,coupling,bond_source,vertex,*edge);
      J = J * fact;
      delta_E += bond_energy_change(moment[boost::source(*edge,graph)],
                                moment[boost::target(*edge,graph)],J,update);
    }
    delta_E += spinfactor[vertex]*site_energy_change(moment[vertex], h, update) * g;
    
    delta_E += onsite_energy_change(moment[vertex], selfinteraction[vertex],update);
    
    // step iv)
    if (delta_E<0 || rng() < exp(-beta*delta_E))
      moment[vertex].update(update); 
  }
}

#endif
