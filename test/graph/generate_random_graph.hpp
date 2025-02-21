/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2015 by Andreas Hehn <hehn@phys.ethz.ch>                          *
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
#ifndef ALPS_TEST_GENERATE_RANDOM_GRAPH_HPP
#define ALPS_TEST_GENERATE_RANDOM_GRAPH_HPP

#include <vector>
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/minmax.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/graph/graph_traits.hpp>

#include <alps/lattice/graphproperties.h>

template <typename Graph, typename RNG>
Graph generate_random_simple_graph(RNG& rng, unsigned int max_vertices, double max_connectivity)
{
    using boost::get;
    boost::random::uniform_int_distribution<unsigned int> rnd_nv(1,max_vertices-1);
    unsigned int num_vertices = rnd_nv(rng)+1;
    unsigned int max_edges = static_cast<unsigned int>(max_connectivity * (num_vertices * num_vertices))+1;
    boost::random::uniform_int_distribution<unsigned int> rnd_edges(0,max_edges-1);
    unsigned int num_edges = rnd_edges(rng)+1;


    unsigned int n = 1;
    std::vector<boost::tuple<unsigned int, unsigned int> > edges;
    for(unsigned int i=0; i < num_edges; ++i)
    {
        unsigned int vtx1 = rng() % n;
        unsigned int vtx2 = vtx1;
        while(vtx2 == vtx1)
            vtx2 = rng() % max_vertices; // vtx2 has to be different from vtx1
        if(vtx2 >= n)
            vtx2 = n++;
        edges.push_back(boost::minmax(vtx1,vtx2));
    }

    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(),edges.end()),edges.end());

    Graph g;
    for(unsigned int i=0; i < edges.size(); ++i)
        add_edge(get<0>(edges[i]), get<1>(edges[i]), g);

    return g;
}

template <typename Graph, typename RNG>
Graph random_isomorphic_graph(RNG& rng, Graph const& g)
{
    using std::swap;
    std::size_t const n = num_vertices(g);
    typename boost::graph_traits<Graph>::vertex_iterator it,end;
    boost::tie(it,end) = vertices(g);
    std::vector<unsigned int> vertex_map(it,end);
    boost::random::uniform_int_distribution<unsigned int> rnd_vtx(0,n-1);
    for(unsigned int i=0; i < n; ++i)
        swap(vertex_map[rnd_vtx(rng)],vertex_map[rnd_vtx(rng)]);

    Graph ng(num_vertices(g));
    typename boost::graph_traits<Graph>::edge_iterator eit,eend;
    for(boost::tie(eit,eend) = edges(g); eit != eend; ++eit)
        add_edge(vertex_map[source(*eit,g)],vertex_map[target(*eit,g)],ng);
    return ng;
}

template <typename Graph, typename RNG>
Graph generate_random_simple_graph_with_colored_edges(RNG& rng, unsigned int max_vertices, double max_connectivity, unsigned int max_colors)
{
    Graph g = generate_random_simple_graph<Graph>(rng, max_vertices, max_connectivity);
    typename boost::property_map<Graph, alps::edge_type_t>::type ecmap = get(alps::edge_type_t(),g);
    boost::random::uniform_int_distribution<unsigned int> rnd_color(0,max_colors-1);
    typename boost::graph_traits<Graph>::edge_iterator it,end;
    for(boost::tie(it,end) = edges(g); it != end; ++it)
        ecmap[*it] = rnd_color(rng);
    return g;
}

#endif // ALPS_TEST_GENERATE_RANDOM_GRAPH_HPP
