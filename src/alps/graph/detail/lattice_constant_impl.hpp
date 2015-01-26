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

#ifndef ALPS_GRAPH_DETAIL_LATTICE_CONSTANT_IMPL_HPP
#define ALPS_GRAPH_DETAIL_LATTICE_CONSTANT_IMPL_HPP

#include <alps/ngs/stacktrace.hpp>

#include <alps/lattice.h>
#include <alps/numeric/vector_functions.hpp>
#include <alps/graph/canonical_properties.hpp>
#include <alps/numeric/matrix.hpp>

#include <boost/array.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>
#include <boost/integer_traits.hpp>

#include <deque>
#include <vector>
#include <cstring>
#include <algorithm>
#include <stdexcept>

#if !defined(USE_COMPRESSED_EMBEDDING) && !defined(USE_COMPRESSED_EMBEDDING2) && !defined(USE_GENERIC_EMBEDDING)
    #define USE_GENERIC_EMBEDDING
#endif

namespace alps {
    namespace graph {

        namespace detail {

            struct embedding_generic_type {

                embedding_generic_type(std::size_t vertices_size, std::size_t edges_size)
                    : hash(0)
                    , counter(new boost::uint8_t(1))
                    , vertices(new std::vector<std::vector<boost::uint16_t> >(vertices_size))
                    , edges(new std::vector<boost::uint64_t>((edges_size >> 6) + ((edges_size & 0x3F) == 0 ? 0 : 1)))
                {}

                embedding_generic_type(embedding_generic_type const & rhs)
                    : hash(rhs.hash)
                    , counter(rhs.counter)
                    , vertices(rhs.vertices)
                    , edges(rhs.edges)
                {
                    assert(*counter < boost::integer_traits<boost::uint8_t>::const_max - 1);
                    ++*counter;
                }

                ~embedding_generic_type() {
                    if (!--*counter) {
                        delete counter;
                        delete vertices;
                        delete edges;
                    }
                }

                bool operator == (embedding_generic_type const & rhs) const {
                    return hash == rhs.hash
                        && *edges == *rhs.edges
                        && *vertices == *rhs.vertices
                    ;
                }

                std::size_t hash;
                boost::uint8_t * counter;
                std::vector<std::vector<boost::uint16_t> > * vertices;
                std::vector<boost::uint64_t> * edges;

                private:
                    embedding_generic_type() {}
            };

            std::size_t hash_value(embedding_generic_type const & value) {
                return value.hash;
            }

            template <typename Graph, typename Lattice> std::vector<std::vector<boost::uint_t<8>::fast> > build_translation_table(
                  Graph const & graph
                , Lattice const & lattice
            ) {
                typedef typename alps::lattice_traits<Lattice>::cell_iterator cell_iterator;
                typedef typename alps::lattice_traits<Lattice>::offset_type offset_type;
                typedef typename alps::lattice_traits<Lattice>::size_type cell_index_type;

                std::vector<std::vector<boost::uint_t<8>::fast> > distance_to_boarder(dimension(lattice), std::vector<boost::uint_t<8>::fast>(num_vertices(graph), num_vertices(graph)));
                std::vector<std::vector<unsigned> > translations(dimension(lattice), std::vector<unsigned>(num_vertices(graph), num_vertices(graph)));
                unsigned vtcs_per_ucell = num_vertices(alps::graph::graph(unit_cell(lattice)));
                for(std::size_t d = 0; d < dimension(lattice); ++d) {
                    for(std::pair<cell_iterator,cell_iterator> c = cells(lattice); c.first != c.second; ++c.first) {
                        offset_type ofst = offset(*c.first,lattice);
                        offset_type move(dimension(lattice));
                        move[d] = -1;
                        std::pair<bool,bool> on_lattice_pbc_crossing = shift(ofst,move,lattice);
                        if(on_lattice_pbc_crossing.first && !on_lattice_pbc_crossing.second) {
                            const cell_index_type cellidx = index(*c.first,lattice);
                            const cell_index_type neighboridx = index(cell(ofst, lattice), lattice);
                            for(unsigned v = 0; v < vtcs_per_ucell; ++v)
                                translations[d][cellidx * vtcs_per_ucell + v] = neighboridx * vtcs_per_ucell + v;
                        }
                    }
                    unsigned v;
                    for (std::vector<unsigned>::const_iterator it = translations[d].begin(); it != translations[d].end(); ++it) {
                        if (*it != num_vertices(graph))
                        {
                            distance_to_boarder[d][v = *it] = 0;
                            while ((v = translations[d][v]) != num_vertices(graph))
                                ++distance_to_boarder[d][*it];
                        }
                    }
                }
                return distance_to_boarder;
            }

            template<typename GeometricInfo, typename Graph, typename Subgraph, typename BreakingVertex> typename boost::disable_if<
                boost::is_same<GeometricInfo, boost::false_type>
            >::type lattice_constant_geometry(
                  GeometricInfo & geometric_info
                , Subgraph const & S
                , Graph const & G
                , std::vector<std::size_t> const & I
                , std::vector<std::vector<boost::uint_t<8>::fast> > const & distance_to_boarder
                , std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const & pinning
                , typename partition_type<Subgraph>::type const & subgraph_orbit
                , BreakingVertex const & breaking_vertex
                , bool inserted
            ) {
                if (inserted)
                    for (typename std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>::const_iterator it = pinning.begin(); it != pinning.end(); ++it) {
                        bool is_valid = true;
                        for(std::size_t d = 0; is_valid && d < distance_to_boarder.size(); ++d)
                            is_valid = distance_to_boarder[d][pinning[breaking_vertex]] <= distance_to_boarder[d][*it];
                        if (is_valid)
                            geometric_info(*it - pinning[breaking_vertex],I[it - pinning.begin()])++;
                    }
            }

            template<typename GeometricInfo, typename Graph, typename Subgraph, typename BreakingVertex> typename boost::enable_if<
                boost::is_same<GeometricInfo, boost::false_type>
            >::type lattice_constant_geometry(
                  GeometricInfo &
                , Subgraph const & S
                , Graph const &
                , std::vector<std::size_t> const &
                , std::vector<std::vector<boost::uint_t<8>::fast> > const &
                , std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const &
                , typename partition_type<Subgraph>::type const &
                , BreakingVertex const &
                , bool
            ) {}

            struct embedding_found {};


            template <typename SubGraph, typename GeometricInfo>
            class lattice_constant_inserter
            {
                typedef SubGraph subgraph_type;
              public:

                template <typename Graph, typename Lattice>
                lattice_constant_inserter(subgraph_type const& S, Graph const& G, Lattice const& L, typename partition_type<subgraph_type>::type const & subgraph_orbit, GeometricInfo& geometric_info, typename boost::graph_traits<subgraph_type>::vertex_descriptor breaking_vertex = 0)
                : subgraph_orbit_(subgraph_orbit), I_(num_vertices(S)), distance_to_boarder_(build_translation_table(G,L)), matches_(), unit_cell_size_(num_vertices(alps::graph::graph(unit_cell(L)))), geometric_info_(geometric_info), breaking_vertex_(breaking_vertex)
                {
                    assert(( boost::is_same<GeometricInfo, boost::false_type>::value ? (get<alps::graph::partition>(canonical_properties(S)) == subgraph_orbit_) : (get<alps::graph::partition>(canonical_properties(S,breaking_vertex_)) == subgraph_orbit_) ));
                    // If the lattice has more than 2 dimensions improve this class
                    assert(distance_to_boarder_.size() < 3);

                    // orbit index => vertices
                    // Io = {(mi, j) : ni element of Vj
                    for (typename partition_type<subgraph_type>::type::const_iterator it = subgraph_orbit_.begin(); it != subgraph_orbit_.end(); ++it)
                        for (typename partition_type<subgraph_type>::type::value_type::const_iterator jt = it->begin(); jt != it->end(); ++jt)
                            I_[*jt] = it - subgraph_orbit_.begin();

                }

                template <typename Graph>
                void operator()(
                      subgraph_type const& S
                    , Graph const& G
                    , std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const & pinning
                ) {
                    assert(( boost::is_same<GeometricInfo, boost::false_type>::value ? (get<alps::graph::partition>(canonical_properties(S)) == subgraph_orbit_) : (get<alps::graph::partition>(canonical_properties(S,breaking_vertex_)) == subgraph_orbit_) ));
                    embedding_generic_type embedding_generic(subgraph_orbit_.size(), num_vertices(S) * (num_vertices(S) + 1) / 2);

                    for (std::vector<std::vector<boost::uint16_t> >::iterator it = embedding_generic.vertices->begin(); it != embedding_generic.vertices->end(); ++it)
                        it->reserve(subgraph_orbit_[it - embedding_generic.vertices->begin()].size());

                    std::size_t bits_per_dim = 0;
                    while ((0x01u << ++bits_per_dim) < num_vertices(S));
                    assert((0x01u << (distance_to_boarder_.size() * bits_per_dim)) < boost::integer_traits<boost::uint16_t>::const_max);

                    std::vector<boost::uint_t<8>::fast> distances(distance_to_boarder_.size(), num_vertices(G));
                    for (typename std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>::const_iterator it = pinning.begin(); it != pinning.end(); ++it)
                        for(std::size_t d = 0; d < distance_to_boarder_.size(); ++d)
                            distances[d] = (std::min)(distances[d], distance_to_boarder_[d][*it]);

                    std::vector<boost::uint16_t> lattice_pinning(pinning.size());
                    for (typename std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>::const_iterator it = pinning.begin(); it != pinning.end(); ++it) {
                        lattice_pinning[it - pinning.begin()] = *it % unit_cell_size_;
                        for(std::size_t d = 0; d < distance_to_boarder_.size(); ++d) {
                            lattice_pinning[it - pinning.begin()] <<= bits_per_dim;
                            lattice_pinning[it - pinning.begin()] += distance_to_boarder_[d][*it] - distances[d];
                        }
                        (*embedding_generic.vertices)[I_[it - pinning.begin()]].push_back(lattice_pinning[it - pinning.begin()]);
                    }
                    for (std::vector<std::vector<boost::uint16_t> >::iterator it = embedding_generic.vertices->begin(); it != embedding_generic.vertices->end(); ++it) {
                        using boost::hash_combine;
                        std::sort(it->begin(), it->end());
                        for (std::vector<boost::uint16_t>::const_iterator jt = it->begin(); jt != it->end(); ++jt)
                            hash_combine(embedding_generic.hash, *jt);
                    }

                    for (std::vector<boost::uint16_t>::iterator it = lattice_pinning.begin(); it != lattice_pinning.end(); ++it) {
                        std::vector<boost::uint16_t>::iterator jt = (*embedding_generic.vertices)[I_[it - lattice_pinning.begin()]].begin();
                        for (; *jt != *it; ++jt);
                        *it = jt - (*embedding_generic.vertices)[I_[it - lattice_pinning.begin()]].begin();
                        for (std::size_t i = 0; i < I_[it - lattice_pinning.begin()]; ++i)
                            *it += (*embedding_generic.vertices)[i].size();
                    }

                    typename boost::graph_traits<subgraph_type>::edge_iterator s_ei, s_ee;
                    for (boost::tie(s_ei, s_ee) = edges(S); s_ei != s_ee; ++s_ei) {
                        std::size_t v1 = std::min(lattice_pinning[source(*s_ei, S)], lattice_pinning[target(*s_ei, S)]);
                        std::size_t v2 = std::max(lattice_pinning[source(*s_ei, S)], lattice_pinning[target(*s_ei, S)]);
                        std::size_t index = v1 * num_vertices(S) - (v1 - 1) * v1 / 2 + v2 - v1;
                        (*embedding_generic.edges)[index >> 6] |= 0x01 << (index & 0x3F);
                    }
                    for (std::vector<boost::uint64_t>::const_iterator it = embedding_generic.edges->begin(); it != embedding_generic.edges->end(); ++it) {
                        using boost::hash_combine;
                        hash_combine(embedding_generic.hash, *it);
                    }

                    lattice_constant_geometry(geometric_info_, S, G, I_, distance_to_boarder_, pinning, subgraph_orbit_, breaking_vertex_, matches_.insert(embedding_generic).second);
                }

                std::size_t get_count() const
                {
                    return matches_.size();
                }
              private:
                typename partition_type<subgraph_type>::type const&             subgraph_orbit_;
                std::vector<std::size_t>                                        I_;
                std::vector<std::vector<boost::uint_t<8>::fast> > const         distance_to_boarder_;
                boost::unordered_set<embedding_generic_type>                    matches_;
                std::size_t const                                               unit_cell_size_;
                GeometricInfo&                                                  geometric_info_;
                typename boost::graph_traits<subgraph_type>::vertex_descriptor const  breaking_vertex_; // ONLY USED WITH GeometricInfo
            };

            struct throw_on_embedding_found
            {
                template <typename SubGraph, typename Graph>
                void operator()(SubGraph const&, Graph const&, std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const&)
                {
                    throw embedding_found();
                }
            };

            template <typename Subgraph>
            struct edge_equal_simple
            {
              private:
                template <typename Graph>
                bool impl(
                  typename boost::graph_traits<Subgraph>::edge_descriptor const & s_e
                , typename boost::graph_traits<Graph>::edge_descriptor const & g_e
                , Subgraph const & S
                , Graph const & G
                , boost::mpl::true_
                ) const {
                    return get(alps::edge_type_t(), S)[s_e] == get(alps::edge_type_t(), G)[g_e];
                }

                template <typename Graph>
                bool impl(
                  typename boost::graph_traits<Subgraph>::edge_descriptor const & s_e
                , typename boost::graph_traits<Graph>::edge_descriptor const & g_e
                , Subgraph const & S
                , Graph const & G
                , boost::mpl::false_
                ) const {
                    return true;
                }

              public:
                template <typename Graph>
                bool operator()(
                  typename boost::graph_traits<Subgraph>::edge_descriptor const & s_e
                , typename boost::graph_traits<Graph>::edge_descriptor const & g_e
                , Subgraph const & S
                , Graph const & G
                ) const {
                    return impl(s_e,g_e,S,G,boost::mpl::bool_<has_property<alps::edge_type_t,Subgraph>::edge_property>());
                }

            };

            template <typename Subgraph>
            struct edge_equal_with_color_symmetries
            {
                typedef typename boost::property_map<Subgraph,alps::edge_type_t>::type::value_type color_type;
                typedef typename color_partition<Subgraph>::type                                color_partition_type;
                typedef std::vector<color_type>                                                 color_map_type;

                static color_type const invalid = boost::integer_traits<color_type>::const_max;

                edge_equal_with_color_symmetries(color_partition_type const& color_partition)
                : color_partition_(color_partition), color_map_(color_partition_.size(),color_type(invalid))
                {
                    BOOST_STATIC_ASSERT(( has_property<alps::edge_type_t,Subgraph>::edge_property ));
                }

                void reset()
                {
                    color_map_.clear();
                    color_map_.resize(color_partition_.size(),invalid);
                }


                template <typename Graph>
                bool operator()(
                      typename boost::graph_traits<Subgraph>::edge_descriptor const & s_e
                    , typename boost::graph_traits<Graph>::edge_descriptor const & g_e
                    , Subgraph const & S
                    , Graph const & G
                ) {
                    color_type const sec = get(alps::edge_type_t(), S)[s_e];
                    color_type const gec = get(alps::edge_type_t(), G)[g_e];
                    assert(sec < color_map_.size());

                    if(color_map_[sec] == invalid && color_partition_[sec] == color_partition_[gec])
                    {
                        // Try to add a mapping from color sec to color gec to the color_map_ while keeping the mapping unique.
                        // (i.e. no two colors sec0,sec1 can be mapped to the same color gec)
                        if(std::find(color_map_.begin(),color_map_.end(), gec) == color_map_.end())
                            color_map_[sec] = gec;
                    }
                    return color_map_[sec] == gec;
                }
            private:
                color_partition_type    color_partition_;
                color_map_type          color_map_;
            };

            template <typename Subgraph>
            struct vertex_equal_simple
            {
              private:
                template <typename Graph>
                bool impl(
                  typename boost::graph_traits<Subgraph>::vertex_descriptor const & s_v
                , typename boost::graph_traits<Graph>::vertex_descriptor const & g_v
                , Subgraph const & S
                , Graph const & G
                , boost::mpl::true_
                ) const {
                    return get(alps::vertex_type_t(), S)[s_v] == get(alps::vertex_type_t(), G)[g_v];
                }

                template <typename Graph>
                bool impl(
                  typename boost::graph_traits<Subgraph>::vertex_descriptor const & s_v
                , typename boost::graph_traits<Graph>::vertex_descriptor const & g_v
                , Subgraph const & S
                , Graph const & G
                , boost::mpl::false_
                ) const {
                    return true;
                }

              public:
                template <typename Graph>
                bool operator()(
                  typename boost::graph_traits<Subgraph>::vertex_descriptor const & s_v
                , typename boost::graph_traits<Graph>::vertex_descriptor const & g_v
                , Subgraph const & S
                , Graph const & G
                ) const {
                    return impl(s_v,g_v,S,G,boost::mpl::bool_<has_property<alps::vertex_type_t,Subgraph>::vertex_property>());
                }
            };

            // TODO: make an object out of walker
            template<typename Subgraph, typename Graph, typename VertexEqual, typename EdgeEqual, typename EmbeddingFoundPolicy> void lattice_constant_walker(
                  typename boost::graph_traits<Subgraph>::vertex_descriptor const & s
                , typename boost::graph_traits<Graph>::vertex_descriptor const & g
                , Subgraph const & S
                , Graph const & G
                , std::deque<std::pair<
                      typename boost::graph_traits<Subgraph>::vertex_descriptor
                    , typename boost::graph_traits<Graph>::vertex_descriptor
                  > > const& queue
                , boost::dynamic_bitset<> const& queued_or_placed
                , boost::dynamic_bitset<> & visited
                , std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> & pinning
                , VertexEqual & vertex_equal
                , EdgeEqual & edge_equal
                , EmbeddingFoundPolicy& register_embedding
            ) {
                typedef typename boost::graph_traits<Subgraph>::vertex_descriptor subgraph_vertex_descriptor;
                typedef typename boost::graph_traits<Graph>::vertex_descriptor       graph_vertex_descriptor;

                // Check if the vertex mapping s->g is valid by checking if...
                // ... the degrees are compatible
                if (out_degree(s, S) > out_degree(g, G))
                    return;
                // ... the vertex types are equal (with respect to symmetries maybe)
                if (!vertex_equal(s, g, S, G))
                    return;
                // ... the existing edges from s are compatible with those of g.
                typename boost::graph_traits<Subgraph>::adjacency_iterator s_ai, s_ae;
                for (boost::tie(s_ai, s_ae) = adjacent_vertices(s, S); s_ai != s_ae; ++s_ai)
                    if (pinning[*s_ai] != num_vertices(G)) {
                        typename boost::graph_traits<Graph>::edge_descriptor e;
                        bool is_e;
                        boost::tie(e, is_e) = edge(g, pinning[*s_ai], G);
                        if (!is_e || !edge_equal( edge(s, *s_ai, S).first , e, S, G) )
                            return;
                    }

                // s->g seems legit => pin s->g.
                visited[g] = true;
                pinning[s] = g;
                // If not all vertices are mapped yet
                if (visited.count() < num_vertices(S)) {
                    // queue mapping adjecent vertices of s
                    std::deque< std::pair<subgraph_vertex_descriptor, graph_vertex_descriptor> > local_queue(queue);
                    boost::dynamic_bitset<> local_queued_or_placed(queued_or_placed);

                    typename boost::graph_traits<Graph>::adjacency_iterator g_ai, g_ae;
                    for (boost::tie(s_ai, s_ae) = adjacent_vertices(s, S); s_ai != s_ae; ++s_ai)
                        if (!local_queued_or_placed[*s_ai]) {
                            local_queued_or_placed[*s_ai] = true;
                            local_queue.push_back(std::make_pair(*s_ai, g));
                        }
                    // take the first entry of the queue
                    // and check if entry.s can be mapped to adjacent vertices of entry.g
                    subgraph_vertex_descriptor t = local_queue[0].first;
                    boost::tie(g_ai, g_ae) = adjacent_vertices(local_queue[0].second, G);
                    local_queue.pop_front();
                    for (; g_ai != g_ae; ++g_ai)
                        if (!visited[*g_ai])
                            detail::lattice_constant_walker(
                                  t
                                , *g_ai
                                , S
                                , G
                                , local_queue
                                , local_queued_or_placed
                                , visited
                                , pinning
                                , vertex_equal
                                , edge_equal
                                , register_embedding
                            );
                } else
                    register_embedding(S,G,pinning);
                pinning[s] = num_vertices(G);
                visited[g] = false;
            }

            // Input: Subgraph, Graph, vertices of G contained in mapping of S on G
            // Output: lattice_constant of S in G containing v
            template<typename Subgraph, typename Graph, typename VertexEqual, typename EdgeEqual, typename EmbeddingFoundPolicy> void lattice_constant_impl(
                  Subgraph const & S
                , Graph const & G
                , std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const & V
                , typename partition_type<Subgraph>::type const & subgraph_orbit
                , VertexEqual & vertex_equal
                , EdgeEqual & edge_equal
                , EmbeddingFoundPolicy& embedding_found_policy
            ) {
                // Assume the vertex desciptor is an unsigned integer type (since we want to use it as an index for a vector)
                BOOST_STATIC_ASSERT((boost::is_unsigned<typename alps::graph_traits<Subgraph>::vertex_descriptor>::value));
                assert(num_vertices(S) > 0);
                // if larger, extend the space
                assert(num_vertices(S) < 21);
                assert(num_edges(S) < 21);

                BOOST_STATIC_ASSERT((boost::is_unsigned<typename alps::graph_traits<Graph>::vertex_descriptor>::value));
                assert(num_vertices(G) > 0);

                // make sure, that a distance in one direction fits in a boost::uint8_t
                assert(std::size_t(num_vertices(G)) < 256 * 256);

                for (typename std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>::const_iterator it = V.begin(); it != V.end(); ++it)
                    for (typename partition_type<Subgraph>::type::const_iterator jt = subgraph_orbit.begin(); jt != subgraph_orbit.end(); ++jt)
                        if (out_degree(jt->front(), S) <= out_degree(*it, G)) {
                            // TODO: shouldn't out_degree be just degree?
                            // TODO: use dynamicbitset
                            boost::dynamic_bitset<> queued_or_placed(num_vertices(S));
                            boost::dynamic_bitset<> visited(num_vertices(G));
                            std::deque<std::pair<
                                  typename boost::graph_traits<Subgraph>::vertex_descriptor
                                , typename boost::graph_traits<Graph>::vertex_descriptor
                            > > queue;
                            std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> pinning(num_vertices(S), num_vertices(G));
                            queued_or_placed[jt->front()] = true;
                            lattice_constant_walker(
                                  jt->front()
                                , *it
                                , S
                                , G
                                , queue
                                , queued_or_placed
                                , visited
                                , pinning
                                , vertex_equal
                                , edge_equal
                                , embedding_found_policy
                            );
                            break;
                        }
            }

        } // end namespace detail
    } // end namespace graph
} // end namespace alps

#endif // ALPS_GRAPH_DETAIL_LATTICE_CONSTANT_IMPL_HPP
