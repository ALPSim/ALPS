/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
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

#ifndef ALPS_GRAPH_LATTICE_CONSTANT
#define ALPS_GRAPH_LATTICE_CONSTANT

#include <alps/ngs/macros.hpp>

#include <alps/lattice/graph_helper.h>
#include <alps/lattice/graphproperties.h>
#include <alps/numeric/vector_functions.hpp>
#include <alps/graph/canonical_properties.hpp>

#include <deque>

namespace alps {
	namespace graph {
	
		namespace detail {

			template<typename Vertex> std::map<unsigned, std::set<Vertex> > lattice_constant_translations(
				  unsigned int invalid
				, std::vector<unsigned int> const & translation
				, std::map<unsigned, std::set<Vertex> > match
			) {
				std::map<unsigned, std::set<Vertex> > moved;
				do {
					moved.clear();
					for (typename std::map<unsigned, std::set<Vertex> >::const_iterator jt = match.begin(); jt != match.end(); ++jt) {
						moved[jt->first] = std::set<Vertex>();
						for (typename std::set<Vertex>::const_iterator kt = jt->second.begin(); kt != jt->second.end(); ++kt)
							if (translation[*kt] == invalid)
								return match;
							else
								moved[jt->first].insert(translation[*kt]);
					}
					swap(match, moved);
				} while(true);
			}

			template<typename Subgraph, typename Graph> void lattice_constant_walker(
				  typename boost::graph_traits<Subgraph>::vertex_descriptor const & s
				, typename boost::graph_traits<Graph>::vertex_descriptor const & g
				, Subgraph const & S
				, Graph const & G
				, std::map<typename boost::graph_traits<Graph>::vertex_descriptor, std::size_t> & I
				, std::set<std::map<unsigned, std::set<typename boost::graph_traits<Graph>::vertex_descriptor> > > & matches
				, std::vector<std::vector<unsigned int> > const & translations
				, std::deque<std::pair<
					  typename boost::graph_traits<Subgraph>::vertex_descriptor
					, typename boost::graph_traits<Graph>::vertex_descriptor
				  > > stack
				, std::map<unsigned, std::set<typename boost::graph_traits<Graph>::vertex_descriptor> > match
				, std::set<typename boost::graph_traits<Subgraph>::vertex_descriptor> placed
				, std::set<typename boost::graph_traits<Graph>::vertex_descriptor> visited
				, std::map<
					  typename boost::graph_traits<Subgraph>::vertex_descriptor
					, typename boost::graph_traits<Graph>::vertex_descriptor
				  > pinning
			) {
				if (out_degree(s, S) > out_degree(g, G))
					return;
				typename boost::graph_traits<Subgraph>::adjacency_iterator s_ai, s_ae;
				for (tie(s_ai, s_ae) = adjacent_vertices(s, S); s_ai != s_ae; ++s_ai)
					if (pinning.find(*s_ai) != pinning.end()) {
						typename boost::graph_traits<Graph>::edge_descriptor e;
						bool is_e;
						tie(e, is_e) = edge(g, pinning[*s_ai], G);
						// TODO: check colored graphs ...
						if (!is_e)
							return;
					}
				// TODO: check colored graphs ...
				visited.insert(g);
				if (match.find(I[s]) == match.end())
					match[I[s]] = std::set<typename boost::graph_traits<Graph>::vertex_descriptor>();
				match[I[s]].insert(g);
				pinning[s] = g;
				if (visited.size() < num_vertices(S)) {
					typename boost::graph_traits<Graph>::adjacency_iterator g_ai, g_ae;
					for (tie(s_ai, s_ae) = adjacent_vertices(s, S); s_ai != s_ae; ++s_ai)
						if (placed.find(*s_ai) == placed.end()) {
							placed.insert(*s_ai);
							stack.push_back(std::make_pair(*s_ai, g));
						}
					typename boost::graph_traits<Subgraph>::vertex_descriptor t = stack[0].first;
					tie(g_ai, g_ae) = adjacent_vertices(stack[0].second, G);
					stack.pop_front();
					for (; g_ai != g_ae; ++g_ai)
						if (visited.find(*g_ai) == visited.end())
							detail::lattice_constant_walker(t, *g_ai, S, G, I, matches, translations, stack, match, placed, visited, pinning);
				} else {
					for (std::vector<std::vector<unsigned int> >::const_iterator it = translations.begin(); it != translations.end(); ++it)
						match = lattice_constant_translations(num_vertices(G), *it, match);
					matches.insert(match);
				}
			}
		}

		// Input: Subgraph, Graph, vertices of G contained in mapping of S on G
		// Output: lattice_constant of S in G containing v
		template<typename Subgraph, typename Graph, typename Lattice, typename LatticeGraph> std::size_t lattice_constant(
			  Subgraph const & S
			, Graph const & G
			, Lattice const & L
			, LatticeGraph const & LG
			, std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> const & v
		) {
			typedef typename alps::graph_helper<Graph>::lattice_type lattice_type;
			typedef typename alps::lattice_traits<lattice_type>::unit_cell_type::graph_type unit_cell_graph_type;

			// TODO: implement that
			if (v.size() != 1)
				ALPS_NGS_THROW_RUNTIME_ERROR("not Impl!")

			// Get the possible translation/coordination vectors from the unit cell
			std::vector<std::vector<int> > coordination_vectors;
			unit_cell_graph_type const unit_cell_graph = alps::graph::graph(unit_cell(L));
			typename boost::graph_traits<unit_cell_graph_type>::edge_iterator it, et;
			for(boost::tie(it, et) = edges(unit_cell_graph); it != et; ++it) {
				using boost::numeric::operators::operator-;
				const std::vector<int> source_offset = boost::get(alps::source_offset_t(), unit_cell_graph, *it);
				const std::vector<int> target_offset = boost::get(alps::target_offset_t(), unit_cell_graph, *it);
				if (target_offset != source_offset)
					coordination_vectors.push_back(target_offset - source_offset);
			}
 
			// Create the translation table if we hit the boundary we set the entry to 'num_vertices' (an invalid vertex index)
			std::vector<std::vector<unsigned int> > translations(
				  coordination_vectors.size()
				, std::vector<unsigned int>(num_vertices(G), num_vertices(G))
			);
			for (boost::tie(it, et) = edges(unit_cell_graph); it != et; ++it) {
				using boost::numeric::operators::operator-;
				// Get cell id (relies on vertex-index = cell-index * vertices_per_cell + vertex-index in cell)
				typename alps::lattice_traits<lattice_type>::offset_type 
					  source_offset = offset(L.cell(source(*it, LG) / num_vertices(unit_cell_graph)), L)
					, target_offset = offset(L.cell(target(*it, LG) / num_vertices(unit_cell_graph)), L)
					, offset_diff = target_offset - source_offset
				;
				for (std::vector<std::vector<int> >::const_iterator jt = coordination_vectors.begin(); jt != coordination_vectors.end(); ++jt) {
					if(offset_diff == *jt)
						translations[jt - coordination_vectors.begin()][source(*it, LG)] = target(*it, LG);
					else if(-offset_diff == *jt)
						translations[jt - coordination_vectors.begin()][target(*it, LG)] = source(*it, LG);
				}
			}

			// orbit index => vertices
			std::set<std::map<unsigned, std::set<typename boost::graph_traits<Graph>::vertex_descriptor> > > matches;
			typename partition_type<Graph>::type orbit = boost::get<2>(canonical_properties(S));
			std::map<typename boost::graph_traits<Subgraph>::vertex_descriptor, std::size_t> I;
			// Io = {(mi, j) : ni element of Vj
			detail::partition_indeces(I, orbit, S);

			for (typename partition_type<Subgraph>::type::const_iterator it = orbit.begin(); it != orbit.end(); ++it)
				if (out_degree(it->front(), S) <= out_degree(v.front(), G)) {
					std::set<typename boost::graph_traits<Subgraph>::vertex_descriptor> placed;
					std::set<typename boost::graph_traits<Graph>::vertex_descriptor> visited;
					std::map<unsigned, std::set<typename boost::graph_traits<Graph>::vertex_descriptor> > match;
					std::deque<std::pair<
						  typename boost::graph_traits<Subgraph>::vertex_descriptor
						, typename boost::graph_traits<Graph>::vertex_descriptor
					> > stack;
					std::map<
						  typename boost::graph_traits<Subgraph>::vertex_descriptor
						, typename boost::graph_traits<Graph>::vertex_descriptor
					> pinning;
					detail::lattice_constant_walker(it->front(), v.front(), S, G, I, matches, translations, stack, match, placed, visited, pinning);
				}
			return matches.size();
		}
	}
}
#endif
