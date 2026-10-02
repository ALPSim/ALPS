/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_GRAPH_LATTICE_CONSTANT_HPP
#define ALPS_GRAPH_LATTICE_CONSTANT_HPP

#include <alps/graph/detail/lattice_constant_impl.hpp>
#include <alps/graph/vertices_of_cell.hpp>

namespace alps {
    namespace graph {

        template<typename Subgraph, typename Graph, typename Lattice> std::size_t lattice_constant(
              Subgraph const & S
            , Graph const & G
            , Lattice const & L
            , typename alps::lattice_traits<Lattice>::cell_descriptor c
        ) {
            typename partition_type<Subgraph>::type subgraph_orbit = boost::get<2>(canonical_properties(S));

            std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> V(vertices_of_cell(c, L));

            detail::vertex_equal_simple<Subgraph> vertex_equal;
            detail::edge_equal_simple<Subgraph>   edge_equal;
            detail::count_translational_invariant_embeddings<Subgraph>  register_embedding(S, G, L, subgraph_orbit);
            detail::lattice_constant_impl(S, G, V, subgraph_orbit, vertex_equal, edge_equal, register_embedding);
            return register_embedding.get_count();
        }

        template<typename Subgraph, typename Graph, typename Lattice> std::size_t lattice_constant(
              Subgraph const & S
            , Graph const & G
            , Lattice const & L
            , typename alps::lattice_traits<Lattice>::cell_descriptor c
            , std::vector<alps::type_type> const& edge_color_mapping
        ) {
            typename partition_type<Subgraph>::type subgraph_orbit = boost::get<2>(canonical_properties(S));

            std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> V(vertices_of_cell(c, L));

            detail::vertex_equal_simple<Subgraph> vertex_equal;
            detail::edge_equal_mapped_colors<Subgraph> edge_equal(edge_color_mapping);
            detail::count_translational_invariant_embeddings<Subgraph>  register_embedding(S, G, L, subgraph_orbit);
            detail::lattice_constant_impl(S, G, V, subgraph_orbit, vertex_equal, edge_equal, register_embedding);
            return register_embedding.get_count();
        }

        template<typename Subgraph, typename Graph, typename Lattice> void lattice_constant(
              alps::numeric::matrix<unsigned int> & lw
            , Subgraph const & S
            , Graph const & G
            , Lattice const & L
            , typename boost::graph_traits<Graph>::vertex_descriptor pin
            , typename boost::graph_traits<Subgraph>::vertex_descriptor b
            , typename partition_type<Subgraph>::type const & subgraph_orbit
        ) {
            assert(get<alps::graph::partition>(canonical_properties(S,b)) == subgraph_orbit);

            detail::vertex_equal_simple<Subgraph> vertex_equal;
            detail::edge_equal_simple<Subgraph>   edge_equal;
            detail::count_how_subgraphs_are_embedded<Subgraph> register_embedding(S, G, L, subgraph_orbit, lw, b);
            detail::lattice_constant_impl_geometric(S, G, pin, subgraph_orbit, vertex_equal, edge_equal, register_embedding, b);
        }
        template<typename Subgraph, typename Graph, typename Lattice> void lattice_constant(
              alps::numeric::matrix<unsigned int> & lw
            , Subgraph const & S
            , Graph const & G
            , Lattice const & L
            , typename boost::graph_traits<Graph>::vertex_descriptor pin
            , typename boost::graph_traits<Subgraph>::vertex_descriptor b
            , typename partition_type<Subgraph>::type const & subgraph_orbit
            , std::vector<alps::type_type> const& edge_color_mapping
        ) {
            assert(get<alps::graph::partition>(canonical_properties(S,b)) == subgraph_orbit);

            detail::vertex_equal_simple<Subgraph> vertex_equal;
            detail::edge_equal_mapped_colors<Subgraph> edge_equal(edge_color_mapping);
            detail::count_how_subgraphs_are_embedded<Subgraph> register_embedding(S, G, L, subgraph_orbit, lw, b);
            detail::lattice_constant_impl_geometric(S, G, pin, subgraph_orbit, vertex_equal, edge_equal, register_embedding, b);
        }
    }
}

#endif //ALPS_GRAPH_LATTICE_CONSTANT_HPP
