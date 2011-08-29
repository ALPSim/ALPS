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

#ifndef ALPS_GRAPH_SUBGRAPH_ITERATOR
#define ALPS_GRAPH_SUBGRAPH_ITERATOR

#include <boost/dynamic_bitset.hpp>
#include <boost/graph/graph_traits.hpp>

#include <set>
#include <deque>

namespace alps {
	namespace graph {
	
		namespace detail {
		
			template<typename Graph> bool is_connected(
				  typename boost::graph_traits<Graph>::vertex_descriptor const & s
				, typename boost::graph_traits<Graph>::vertex_descriptor const & t
				, Graph const & G
			) {
				std::set<typename boost::graph_traits<Graph>::vertex_descriptor> V;
				std::deque<typename boost::graph_traits<Graph>::vertex_descriptor> S(1, s);
				typename boost::graph_traits<Graph>::adjacency_iterator ai, ae;
				while (S.size()) {
					for (tie(ai, ae) = adjacent_vertices(S.front(), G); ai != ae; ++ai)
						if (*ai == t)
							return true;
						else if (V.insert(*ai).second)
							S.push_back(*ai);
					S.pop_front();
				}
				return false;
			}

			template<typename Graph> void subgraphs_helper(
				  std::set<boost::dynamic_bitset<> > & L
				, typename boost::graph_traits<Graph>::vertex_descriptor const & s
				, typename boost::graph_traits<Graph>::vertex_descriptor const & t
				, Graph G
			) {
				remove_edge(s, t, G);
				if (num_edges(G)) {
					if (!out_degree(s, G) || !out_degree(t, G) || is_connected(s, t, G)) {
						boost::dynamic_bitset<> l(num_vertices(G) * (num_vertices(G) + 1) / 2);
						typename boost::graph_traits<Graph>::edge_iterator ei, ee;
						for	(boost::tie(ei, ee) = edges(G); ei != ee; ++ei)
							l[source(*ei, G) * num_vertices(G) - (source(*ei, G) - 1) * source(*ei, G) / 2 + target(*ei, G) - source(*ei, G)] = true;
						if (L.insert(l).second) {
							typename boost::graph_traits<Graph>::edge_iterator ei, ee;
							for	(boost::tie(ei, ee) = edges(G); ei != ee; ++ei)
								detail::subgraphs_helper(L, source(*ei, G), target(*ei, G), G);
						}
					}
				}
			}

		}

		template<typename Graph> void subgraphs(std::set<boost::dynamic_bitset<> > & L, Graph const & G) {
			typename boost::graph_traits<Graph>::edge_iterator it, end;
			for (boost::tie(it, end) = edges(G); it != end; ++it)
				detail::subgraphs_helper(L, source(*it, G), target(*it, G), G);
		}

	}
}

#endif
