/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2015 by Lukas Gamper <gamperl@gmail.com>                   *
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

#ifndef ALPS_GRAPH_CANONICAL_PROPERTIES_HPP
#define ALPS_GRAPH_CANONICAL_PROPERTIES_HPP

#include <alps/graph/canonical_properties_traits.hpp>
#include <alps/graph/detail/canonical_properties_impl.hpp>
#include <boost/graph/properties.hpp>

#include <boost/mpl/if.hpp>
#include <boost/container/flat_map.hpp>

#ifndef NDEBUG
#include <alps/graph/detail/assert_helpers.hpp>
#endif //NDEBUG

namespace alps {
    namespace graph {
        // McKay’s canonical isomorph function Cm(G) is deﬁned to be
        // Cm(G) = max{ Gpi: (pi, nu) is a leaf of T(G) }
        // Input: graph G
        // Output: canonical ordering, canonical label Gpi and orbits of the vertices under the action of the automorphism group of G
        template<typename Graph>
        typename canonical_properties_type<Graph>::type
        canonical_properties(Graph const & G) {
            // The McKay Algorithm works for simple graphs only!
            // (Note: Edges connected to the same vertex on both sides might work,
            // but were not considered in the papers and not tested here.)
            assert( detail::assert_helpers::is_simple_graph(G) );
            typename partition_type<Graph>::type pi;
            // pi = (V1, V2, ..., Vr), Vi = (n1, n2, ..., nk), ni element of G
            // uncolored graphs: pi is a unit partition (pi has only one part)
            // vertex colored graphs: pi has for each color one part
            detail::initial_partition(G, pi, boost::mpl::bool_<has_property<alps::vertex_type_t, Graph>::vertex_property>());
            typedef typename boost::mpl::if_c<
                  has_property<alps::vertex_type_t, Graph>::vertex_property
                , detail::label::simple_vertex_coloring_policy
                , detail::label::no_coloring_policy
            >::type vertex_policy;
            typedef typename boost::mpl::if_c<
                  has_property<alps::edge_type_t, Graph>::edge_property
                , detail::label::simple_edge_coloring_policy
                , detail::label::no_coloring_policy
            >::type edge_policy;

            detail::label::graph_label_creator<Graph,vertex_policy,edge_policy> label_creator(G);
            // create canonical properties
            return detail::canonical_properties_impl(G, pi, label_creator);
        }

        // McKay’s canonical isomorph function Cm(G) is deﬁned to be
        // Cm(G) = max{ Gpi: (pi, nu) is a leaf of T(G) }
        // Input: graph G
        // Output: canonical ordering, canonical label Gpi and orbits of the vertices under the action of the automorphism group of G
        // Optionally outputs the map mapping (color_mapping) the colors of the input graph to
        // its canonical representative with respect to the edge color symmetries.
        template<typename Graph>
        typename canonical_properties_type<Graph>::type
        canonical_properties(Graph const & G, typename color_partition<Graph>::type const& c, boost::container::flat_map<alps::type_type, alps::type_type> * const color_mapping = NULL) {
            // The McKay Algorithm works for simple graphs only!
            // (Note: Edges connected to the same vertex on both sides might work,
            // but were not considered in the papers and not tested here.)
            assert( detail::assert_helpers::is_simple_graph(G) );
            typename partition_type<Graph>::type pi;
            // pi = (V1, V2, ..., Vr), Vi = (n1, n2, ..., nk), ni element of G
            // uncolored graphs: pi is a unit partition (pi has only one part)
            // vertex colored graphs: pi has for each color one part
            detail::initial_partition(G, pi, boost::mpl::bool_<has_property<alps::vertex_type_t, Graph>::vertex_property>());

            // Build the label coloring helpers
            BOOST_STATIC_ASSERT(( has_property<alps::edge_type_t, Graph>::edge_property ));
            typedef typename boost::mpl::if_c<
                  has_property<alps::vertex_type_t, Graph>::vertex_property
                , detail::label::simple_vertex_coloring_policy
                , detail::label::no_coloring_policy
            >::type vertex_policy;
            detail::label::graph_label_creator<Graph, vertex_policy, detail::label::edge_coloring_with_symmetries_policy> label_creator(G);
            label_creator.set_color_partition(c);
            // create canonical properties
            typename canonical_properties_type<Graph>::type r(detail::canonical_properties_impl(G, pi, label_creator));
            if(color_mapping != NULL)
                *color_mapping = label_creator.get_color_mapping();
            return r;
        }

        // Function overload for the previous function to ease generic programming.
        template<typename Graph>
        typename canonical_properties_type<Graph>::type
        canonical_properties(Graph const & G, no_color_symmetry) {
            return canonical_properties(G);
        }

        // McKay’s canonical isomorph function Cm(G) is deﬁned to be
        // Cm(G) = max{ Gpi: (pi, nu) is a leaf of T(G) }
        // Input: graph G, symmetry breaking vertex v
        // Output: canonical ordering, canonical label Gpi and orbits of the vertices under the action of the automorphism group of G
        template<typename Graph>
        typename canonical_properties_type<Graph>::type
        canonical_properties(Graph const & G, typename boost::graph_traits<Graph>::vertex_descriptor v) {
            // The McKay Algorithm works for simple graphs only!
            // (Note: Edges connected to the same vertex on both sides might work,
            // but were not considered in the papers and not tested here.)
            assert( detail::assert_helpers::is_simple_graph(G) );
            assert( detail::assert_helpers::graph_has_vertex(G,v));
            typename partition_type<Graph>::type pi;
            // pi = (V1, V2, ..., Vr), Vi = (n1, n2, ..., nk), ni element of G
            // uncolored graphs: pi is a unit partition (pi has only one part)
            // vertex colored graphs: pi has for each color one part
            if(num_vertices(G) == 1)
                detail::initial_partition(G, pi, boost::mpl::bool_<has_property<alps::vertex_type_t, Graph>::vertex_property>());
                // symmetry breaking vertex is the only vertex -> there is only a single partition
            else
                detail::initial_partition(G, pi, v, boost::mpl::bool_<has_property<alps::vertex_type_t, Graph>::vertex_property>());
            typedef typename boost::mpl::if_c<
                  has_property<alps::vertex_type_t, Graph>::vertex_property
                , detail::label::simple_vertex_coloring_policy
                , detail::label::no_coloring_policy
            >::type vertex_policy;
            typedef typename boost::mpl::if_c<
                  has_property<alps::edge_type_t, Graph>::edge_property
                , detail::label::simple_edge_coloring_policy
                , detail::label::no_coloring_policy
            >::type edge_policy;

            detail::label::graph_label_creator<Graph,vertex_policy,edge_policy> label_creator(G);
            // create canonical properties
            return detail::canonical_properties_impl(G, pi, label_creator);
        }
    } // end namespace graph
} // end namespace alps

#endif // ALPS_GRAPH_CANONICAL_PROPERTIES_HPP
