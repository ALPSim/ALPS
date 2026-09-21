/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013 by Andreas Hehn <hehn@phys.ethz.ch>                          *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "generate_random_graph.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <alps/lattice.h>
#include <alps/graph/canonical_properties.hpp>
#include <boost/random/mersenne_twister.hpp>

template <typename RNG, typename Graph>
bool canonical_label_test(RNG& rng, Graph const& g, unsigned int iterations = 100)
{
    using alps::graph::canonical_properties;
    typedef typename alps::graph::canonical_properties_type<Graph>::type canonical_properties_type;

    canonical_properties_type const gp(canonical_properties(g));
    bool ok = true;
    for(unsigned int i = 0; i < iterations; ++i)
    {
       Graph const ng = random_isomorphic_graph(rng,g);
       canonical_properties_type const ngp(canonical_properties(ng));
       ok = ok && (get<alps::graph::label>(gp) == get<alps::graph::label>(ngp));
       if(!ok)
           std::cout << "ERROR! graph label mismatch: "
               << get<alps::graph::label>(gp)
               << " != " << get<alps::graph::label>(ngp)
               << std::endl;
    }
    return ok;
}

int main()
{
    boost::random::mt19937 rng(23);

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_type;

    bool ok = true;
    for(unsigned int i=0; i < 50; ++i)
    {
        graph_type g = generate_random_simple_graph<graph_type>(rng, 20, 0.6);
        std::cout << "#v = " << num_vertices(g) << " #e = " << num_edges(g) << "\t";
        std::cout << get<alps::graph::label>(alps::graph::canonical_properties(g)) << std::endl;
        ok = ok && canonical_label_test(rng,g);
    }
    return ok ? 0 : -1;
}
