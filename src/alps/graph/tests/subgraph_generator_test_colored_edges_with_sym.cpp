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

#include <gtest/gtest.h>
#include <numeric>
#include "graph_assertions.hpp"
#include <alps/parser/xslt_path.h>
#include <alps/graph/subgraph_generator.hpp>
#include <alps/graph/utils.hpp>
#include <boost/graph/adjacency_list.hpp>

static unsigned int const test_graph_size = 7;

template <typename Graph>
void subgraph_generator_with_color_symmetries_test(unsigned int order)
{
    std::ifstream in(alps::search_xml_library_path("lattices.xml"));
    alps::Parameters parm;
    parm["LATTICE"] = "anisotropic triangular lattice";
    parm["L"]       = 2*order+1;
    alps::graph_helper<> alps_lattice(in,parm);

    typedef alps::coordinate_graph_type lattice_graph_type;
    lattice_graph_type& lattice_graph = alps_lattice.graph();
    std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> pin(1, 2*order*order);

    typename alps::graph::color_partition<Graph>::type color_sym_group;
    color_sym_group[0] = 0;
    color_sym_group[1] = 0;
    color_sym_group[2] = 0;
    std::vector<std::pair<Graph, typename alps::graph::canonical_properties_type<Graph>::type> > v = alps::graph::generate_subgraphs(Graph(), lattice_graph, pin, order, color_sym_group);
    std::vector<std::size_t> counts(order + 1, 0);
    for (auto current = v.begin(); current != v.end(); ++current) {
        const auto edges = num_edges(current->first);
        ASSERT_LE(edges, order);
        ++counts[edges];
    }
    EXPECT_EQ(counts, (std::vector<std::size_t>{1, 1, 2, 7, 23, 104, 539, 3055}));
    EXPECT_EQ(std::accumulate(counts.begin(), counts.end(), std::size_t(0)), 3732u);
}

TEST(GraphSubgraphs, subgraph_generator_test_colored_edges_with_sym)
{
    typedef boost::adjacency_list<boost::vecS, boost::vecS,boost::undirectedS, boost::no_property, boost::property<alps::edge_type_t,alps::type_type> > graph_type;
    subgraph_generator_with_color_symmetries_test<graph_type>(test_graph_size);

}
