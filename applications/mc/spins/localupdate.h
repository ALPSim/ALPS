/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1999-2006 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabian Stoeckli <fabstoec@student.ethz.ch>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
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

template <class Graph, class MomentMap, class RNG, class CouplingMap, class SpinfactorMap>
void local_update(const Graph& graph, MomentMap& moment,
                  double beta, RNG& rng,  const CouplingMap& coupling,
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
      MAT J = coupling[*edge];
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
      MAT J = coupling[*edge];
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
