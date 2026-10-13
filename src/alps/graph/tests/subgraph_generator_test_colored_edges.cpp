/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                              Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//#define USE_COMPRESSED_EMBEDDING2

#include <gtest/gtest.h>
#include <numeric>
#include "graph_assertions.hpp"
#include <alps/parser/xslt_path.h>
#include <alps/graph/subgraph_generator.hpp>
#include <alps/graph/utils.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <iostream>



enum { test_graph_size = 7 };

template <typename Graph>
void subgraph_generator_test(unsigned int order_ )
{
    std::ifstream in(alps::search_xml_library_path("lattices.xml"));
    alps::Parameters parm;
    parm["LATTICE"] = "coupled ladders";
    parm["L"] = 2*order_+1;
    alps::graph_helper<> alps_lattice(in,parm);

    typedef alps::coordinate_graph_type lattice_graph_type;
    lattice_graph_type& lattice_graph = alps_lattice.graph();
    std::map<unsigned int, unsigned int> edge_type_map;
    edge_type_map[0] = 0;
    edge_type_map[1] = 1;
    edge_type_map[2] = 0;
    alps::graph::remap_edge_types(lattice_graph,edge_type_map);


    typedef alps::graph::subgraph_generator<Graph,lattice_graph_type> graph_gen_type;
    std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> pin(1, 2*order_*order_);
    graph_gen_type graph_gen(lattice_graph, pin);

    typename graph_gen_type::iterator it,end;
    boost::tie(it,end) = graph_gen.generate_up_to_n_edges(order_);
    std::vector<std::size_t> counts(order_ + 1, 0);
    for (auto current = it; current != end; ++current) {
        const auto edges = num_edges(current->first);
        ASSERT_LE(edges, order_);
        ++counts[edges];
    }
    EXPECT_EQ(counts, (std::vector<std::size_t>{1, 2, 2, 6, 13, 35, 100, 309}));
    EXPECT_EQ(std::accumulate(counts.begin(), counts.end(), std::size_t(0)), 468u);
}

TEST(GraphSubgraphs, subgraph_generator_test_colored_edges)
{
    typedef boost::adjacency_list<boost::vecS, boost::vecS,boost::undirectedS, boost::no_property, boost::property<alps::edge_type_t,alps::type_type> > graph_type;
    subgraph_generator_test<graph_type>(test_graph_size);

}
