/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
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

#ifndef ALPS_GRAPH_SUBGRAPHS
#define ALPS_GRAPH_SUBGRAPHS

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
                    for (boost::tie(ai, ae) = adjacent_vertices(S.front(), G); ai != ae; ++ai)
                        if (*ai == t)
                            return true;
                        else if (V.insert(*ai).second)
                            S.push_back(*ai);
                    S.pop_front();
                }
                return false;
            }
            
            template<typename Graph> bool no_disconnect(
                  typename boost::graph_traits<Graph>::vertex_descriptor const & vd
                , Graph G
            ) {
                typename boost::graph_traits<Graph>::adjacency_iterator ai, aj, ae;
                boost::tie(ai, ae) = adjacent_vertices(vd, G);
                clear_vertex( vd, G );
                for( aj = ai++ ; ai != ae; ++ai )
                {
                    if( !is_connected( *ai, *aj, G ) )
                        return false;
                    ++aj;
                }
                return true;
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
                        for (boost::tie(ei, ee) = edges(G); ei != ee; ++ei) {
                            typename boost::graph_traits<Graph>::vertex_descriptor v1 = std::min(source(*ei,G),target(*ei,G));
                            typename boost::graph_traits<Graph>::vertex_descriptor v2 = std::max(source(*ei,G),target(*ei,G));
                            l[v1 * num_vertices(G) - (v1 - 1) * v1 / 2 + v2 - v1] = true;
                        }
                        if (L.insert(l).second) {
                            typename boost::graph_traits<Graph>::edge_iterator ei, ee;
                            for(boost::tie(ei, ee) = edges(G); ei != ee; ++ei)
                                detail::subgraphs_helper(L, source(*ei, G), target(*ei, G), G);
                        }
                    }
                }
            }
            
            template<typename Graph> void subgraphs_helper_strong(
                  std::set<boost::dynamic_bitset<> > & L
                , typename boost::graph_traits<Graph>::vertex_descriptor const & vd
                , Graph G
            ) {
                if (num_edges(G)>1)
                {
                    if ( out_degree(vd, G)==1 || no_disconnect(vd, G) )
                    {
                       clear_vertex(vd, G);
                       boost::dynamic_bitset<> l(num_vertices(G) * (num_vertices(G) + 1) / 2);
                       typename boost::graph_traits<Graph>::edge_iterator ei, ee;
                       for (boost::tie(ei, ee) = edges(G); ei != ee; ++ei) 
                       {
                           typename boost::graph_traits<Graph>::vertex_descriptor v1 = std::min(source(*ei,G),target(*ei,G));
                           typename boost::graph_traits<Graph>::vertex_descriptor v2 = std::max(source(*ei,G),target(*ei,G));
                           l[v1 * num_vertices(G) - (v1 - 1) * v1 / 2 + v2 - v1] = true;
                       }
                       if (L.insert(l).second) 
                       {
                            typename boost::graph_traits<Graph>::vertex_iterator vi, ve;
                            for(boost::tie(vi, ve) = vertices(G); vi != ve; ++vi)
                                if( out_degree(*vi, G) )
                                    detail::subgraphs_helper_strong(L, *vi, G);
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
        
        template<typename Graph> void subgraphs_strong(std::set<boost::dynamic_bitset<> > & L, Graph const & G) {
            typename boost::graph_traits<Graph>::vertex_iterator it, end;
            for (boost::tie(it, end) = vertices(G); it != end; ++it)
                detail::subgraphs_helper_strong(L, *it, G);
        }

    }
}

#endif
