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

#include <alps/graph/subgraph_generator.hpp>
#include <alps/graph/utils.hpp>
#include <boost/graph/adjacency_list.hpp>

static unsigned int const test_graph_size = 7;

template <typename Graph>
void subgraph_generator_with_color_symmetries_test(unsigned int order)
{
    std::ifstream in("../../lib/xml/lattices.xml");
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
    std::cout << v.size() << std::endl;
}

int main()
{
    typedef boost::adjacency_list<boost::vecS, boost::vecS,boost::undirectedS, boost::no_property, boost::property<alps::edge_type_t,alps::type_type> > graph_type;
    subgraph_generator_with_color_symmetries_test<graph_type>(test_graph_size);
    return 0;
}
