/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef ALPS_GRAPH_IS_EMBEDDABLE_HPP
#define ALPS_GRAPH_IS_EMBEDDABLE_HPP

#ifdef USE_LATTICE_CONSTANT_2D
#include <alps/graph/lattice_constant_2d.hpp>
#define ALPS_GRAPH_LATTICE_CONSTANT_HPP
#define ALPS_GRAPH_IS_EMBEDDABLE_HPP
#endif // USE_LATTICE_CONSTANT_2D

#include <alps/graph/detail/lattice_constant_impl.hpp>
#include <alps/graph/utils.hpp>

namespace alps {
    namespace graph {

        template<typename Subgraph, typename Graph> bool is_embeddable(
              Subgraph const & S
            , Graph const & G
            , std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const& V
            , typename partition_type<Subgraph>::type const & subgraph_orbit
        ) {
            assert( detail::assert_helpers::graph_has_vertices(G,V) );
            assert(get<alps::graph::partition>(canonical_properties(S)) == subgraph_orbit);
            std::vector<std::vector<boost::uint_t<8>::fast> > distance_to_boarder;

            detail::vertex_equal_simple<Subgraph>   vertex_equal;
            detail::edge_equal_simple<Subgraph>     edge_equal;
            detail::throw_on_embedding_found        throw_on_found_embedding;
            try {
                detail::lattice_constant_impl(S, G, V, subgraph_orbit, vertex_equal, edge_equal, throw_on_found_embedding);
                return false;
            } catch (detail::embedding_found e) {
                return true;
            }
        }

        template<typename Subgraph, typename Graph> bool is_embeddable(
              Subgraph const & S
            , Graph const & G
            , std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const& V
            , typename partition_type<Subgraph>::type const & subgraph_orbit
            , typename color_partition<Subgraph>::type const & color_partition
        ) {
            assert( false || "FIXME" );
            assert( detail::assert_helpers::graph_has_vertices(G,V) );
            assert(get<alps::graph::partition>(canonical_properties(S)) == subgraph_orbit);
            std::vector<std::vector<boost::uint_t<8>::fast> > distance_to_boarder;

            detail::vertex_equal_simple<Subgraph>   vertex_equal;
            detail::throw_on_embedding_found        throw_on_found_embedding;
            std::vector< std::vector<alps::type_type> > color_mappings = get_all_color_mappings_from_color_partition(S, color_partition);
            for(std::vector< std::vector<alps::type_type> >::iterator it = color_mappings.begin(); it != color_mappings.end(); ++it)
            { 
                try {
                    detail::edge_equal_mapped_colors<Subgraph> edge_equal(*it);
                    detail::lattice_constant_impl(S, G, V, subgraph_orbit, vertex_equal, edge_equal, throw_on_found_embedding);
                } catch (detail::embedding_found e) {
                    return true;
                }
            }
            return false;
        }

        template<typename Subgraph, typename Graph> bool is_embeddable(
              Subgraph const & S
            , Graph const & G
            , typename partition_type<Subgraph>::type const & subgraph_orbit
        ) {
            assert(get<alps::graph::partition>(canonical_properties(S)) == subgraph_orbit);
            std::vector<std::vector<boost::uint_t<8>::fast> > distance_to_boarder;

            detail::vertex_equal_simple<Subgraph>   vertex_equal;
            detail::edge_equal_simple<Subgraph>     edge_equal;
            detail::throw_on_embedding_found                    throw_on_found_embedding;
            try {
                typename boost::graph_traits<Graph>::vertex_iterator vt, ve;
                std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> V;
                for (boost::tie(vt, ve) = vertices(G); vt != ve; ++vt) {
                    V.clear();
                    V.push_back(*vt);
                    detail::lattice_constant_impl(S, G, V, subgraph_orbit, vertex_equal, edge_equal, throw_on_found_embedding);
                }
                return false;
            } catch (detail::embedding_found e) {
                return true;
            }
        }

        /**
          * alps::edge_type_t must be an integer value which can be used as index in a vector.
          */
        template<typename Subgraph, typename Graph> bool is_embeddable(
              Subgraph const & S
            , Graph const & G
            , typename partition_type<Subgraph>::type const & subgraph_orbit
            , typename color_partition<Subgraph>::type const & color_partition
        ) {
            assert(get<alps::graph::partition>(canonical_properties(S)) == subgraph_orbit);
            std::vector<std::vector<boost::uint_t<8>::fast> > distance_to_boarder;

            detail::vertex_equal_simple<Subgraph>               vertex_equal;
            detail::throw_on_embedding_found                    throw_on_found_embedding;

            std::vector< std::vector<alps::type_type> > color_mappings = get_all_color_mappings_from_color_partition(S, color_partition);

            try {
                typename boost::graph_traits<Graph>::vertex_iterator vt, ve;
                std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> V;
                for (boost::tie(vt, ve) = vertices(G); vt != ve; ++vt) {
                    V.clear();
                    V.push_back(*vt);
                    for(std::vector< std::vector<alps::type_type> >::iterator it = color_mappings.begin(); it != color_mappings.end(); ++it)
                    { 
                        detail::edge_equal_mapped_colors<Subgraph> edge_equal(*it);
                        detail::lattice_constant_impl(S, G, V, subgraph_orbit, vertex_equal, edge_equal, throw_on_found_embedding);
                    }
                }
            } catch (detail::embedding_found e) {
                return true;
            }
            return false;
        }

    }
}

#endif // ALPS_GRAPH_IS_EMBEDDABLE_HPP
