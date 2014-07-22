/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifdef USE_LATTICE_CONSTANT_2D
#include <alps/graph/lattice_constant_2d.hpp>
#define ALPS_GRAPH_LATTICE_CONSTANT_HPP
#define ALPS_GRAPH_IS_EMBEDDABLE_HPP
#endif // USE_LATTICE_CONSTANT_2D

#ifndef ALPS_GRAPH_LATTICE_CONSTANT_HPP
#define ALPS_GRAPH_LATTICE_CONSTANT_HPP

#include <alps/graph/detail/lattice_constant_impl.hpp>

namespace alps {
    namespace graph {

        template<typename Subgraph, typename Graph, typename Lattice> std::size_t lattice_constant(
              Subgraph const & S
            , Graph const & G
            , Lattice const & L
            , typename alps::lattice_traits<Lattice>::cell_descriptor c
        ) {
            typename partition_type<Subgraph>::type subgraph_orbit = boost::get<2>(canonical_properties(S));

            typename alps::lattice_traits<Lattice>::size_type const cell_id = index(c, L);
            std::size_t unit_cell_size = num_vertices(alps::graph::graph(unit_cell(L)));
            std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> V;
            for(unsigned v = 0; v < unit_cell_size; ++v)
                V.push_back(cell_id * unit_cell_size + v);

            detail::vertex_equal_simple<Subgraph> vertex_equal;
            detail::edge_equal_simple<Subgraph>   edge_equal;
            boost::false_type no_argument;
            detail::lattice_constant_inserter<Subgraph, boost::false_type>    register_embedding(S, G, L, subgraph_orbit, no_argument);
            detail::lattice_constant_impl(S, G, V, subgraph_orbit, vertex_equal, edge_equal, register_embedding);
            return register_embedding.get_count();
        }

        template<typename Subgraph, typename Graph, typename Lattice> std::size_t lattice_constant(
              alps::numeric::matrix<unsigned int> & lw
            , Subgraph const & S
            , Graph const & G
            , Lattice const & L
            , typename alps::lattice_traits<Lattice>::cell_descriptor c
            , typename boost::graph_traits<Subgraph>::vertex_descriptor b
            , typename partition_type<Subgraph>::type const & subgraph_orbit
        ) {
            assert(get<alps::graph::partition>(canonical_properties(S,b)) == subgraph_orbit);

            typename alps::lattice_traits<Lattice>::size_type const cell_id = index(c, L);
            std::size_t unit_cell_size = num_vertices(alps::graph::graph(unit_cell(L)));
            std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> V;
            for(unsigned v = 0; v < unit_cell_size; ++v)
                V.push_back(cell_id * unit_cell_size + v);

            detail::vertex_equal_simple<Subgraph> vertex_equal;
            detail::edge_equal_simple<Subgraph>   edge_equal;
            detail::lattice_constant_inserter<Subgraph, alps::numeric::matrix<unsigned int> >    register_embedding(S, G, L, subgraph_orbit, lw, b);
            detail::lattice_constant_impl(S, G, V, subgraph_orbit, vertex_equal, edge_equal, register_embedding);
            return register_embedding.get_count();
        }
    }
}

#endif //ALPS_GRAPH_LATTICE_CONSTANT_HPP
