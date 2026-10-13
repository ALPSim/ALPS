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
#include "graph_assertions.hpp"
#include <alps/lattice.h> // for cout << graph
#include "generate_random_graph.hpp"
#include <alps/graph/canonical_properties.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/algorithm/cxx11/is_sorted.hpp>


using alps::graph::canonical_properties;

typedef boost::property<alps::edge_type_t,unsigned int> edge_props;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, edge_props> graph_type;
typedef boost::graph_traits<graph_type>::edge_descriptor edge_descriptor;
typedef boost::property_map<graph_type,alps::edge_type_t>::type edge_color_map_type;



TEST(GraphOrbits, orbit_test1) {
    // 0---3
    // |   *
    // 1***2
    graph_type g(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(1, 2, g).first;
        edge_color[e] = 1;
        e = add_edge(2, 3, g).first;
        edge_color[e] = 1;
        e = add_edge(0, 3, g).first;
        edge_color[e] = 0;
    }

    alps::graph::canonical_properties_type<graph_type>::type gp(canonical_properties(g));

    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0;
    color_symmetry[1] = 0;

    alps::graph::canonical_properties_type<graph_type>::type gp_with_sym(canonical_properties(g,color_symmetry));

    expect_partition(get<alps::graph::partition>(gp), {{0}, {1, 3}, {2}});
    expect_partition(get<alps::graph::partition>(gp_with_sym), {{0}, {1, 3}, {2}});
}

TEST(GraphOrbits, orbit_test2) {
    //    2
    //  /  .
    // 0++++1
    //  .  /
    //   3
    graph_type g(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(0, 2, g).first;
        edge_color[e] = 1;
        e = add_edge(0, 3, g).first;
        edge_color[e] = 2;
        e = add_edge(1, 2, g).first;
        edge_color[e] = 2;
        e = add_edge(1, 3, g).first;
        edge_color[e] = 1;
    }

    alps::graph::canonical_properties_type<graph_type>::type gp(canonical_properties(g));

    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0;
    color_symmetry[1] = 0;
    color_symmetry[2] = 0;

    alps::graph::canonical_properties_type<graph_type>::type gp_with_sym(canonical_properties(g,color_symmetry));

    expect_partition(get<alps::graph::partition>(gp), {{0, 1}, {2, 3}});
    expect_partition(get<alps::graph::partition>(gp_with_sym), {{0, 1}, {2, 3}});
}

TEST(GraphOrbits, orbit_test3) {
    //    2+++4
    //  /    /
    // 0++++1
    //  .
    //   3
    graph_type g(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(0, 2, g).first;
        edge_color[e] = 1;
        e = add_edge(0, 3, g).first;
        edge_color[e] = 2;
        e = add_edge(1, 4, g).first;
        edge_color[e] = 1;
        e = add_edge(2, 4, g).first;
        edge_color[e] = 0;
    }

    alps::graph::canonical_properties_type<graph_type>::type gp(canonical_properties(g));

    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0;
    color_symmetry[1] = 0;
    color_symmetry[2] = 0;

    alps::graph::canonical_properties_type<graph_type>::type gp_with_sym(canonical_properties(g,color_symmetry));

    expect_partition(get<alps::graph::partition>(gp), {{0}, {1}, {2}, {3}, {4}});
    expect_partition(get<alps::graph::partition>(gp_with_sym), {{0}, {1}, {2}, {3}, {4}});
}

TEST(GraphOrbits, orbit_test4) {
    {
        //     2
        //     |
        // 4---0+++3
        //     |
        //     1

        graph_type g(4);
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(0, 2, g).first;
        edge_color[e] = 0;
        e = add_edge(0, 3, g).first;
        edge_color[e] = 1;
        e = add_edge(0, 4, g).first;
        edge_color[e] = 0;

        alps::graph::canonical_properties_type<graph_type>::type gp(canonical_properties(g));
        expect_partition(get<alps::graph::partition>(gp), {{0}, {1, 2, 4}, {3}});
    }

    {
        //     2
        //     |
        // 4+++0+++3
        //     |
        //     1
        graph_type g(4);
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(0, 2, g).first;
        edge_color[e] = 0;
        e = add_edge(0, 3, g).first;
        edge_color[e] = 1;
        e = add_edge(0, 4, g).first;
        edge_color[e] = 1;

        alps::graph::canonical_properties_type<graph_type>::type gp(canonical_properties(g));
        expect_partition(get<alps::graph::partition>(gp), {{0}, {1, 2}, {3, 4}});
    }
}

TEST(GraphOrbits, orbit_test5) {
    //
    //  7---0---6---2       // c0 ---
    //  +           +       // c1 +++
    //  1---5---3---4
    //
    {
        graph_type g(8);
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 7, g).first;
        edge_color[e] = 0;
        e = add_edge(0, 6, g).first;
        edge_color[e] = 0;
        e = add_edge(2, 6, g).first;
        edge_color[e] = 0;
        e = add_edge(2, 4, g).first;
        edge_color[e] = 1;
        e = add_edge(3, 4, g).first;
        edge_color[e] = 0;
        e = add_edge(3, 5, g).first;
        edge_color[e] = 0;
        e = add_edge(1, 5, g).first;
        edge_color[e] = 0;
        e = add_edge(1, 7, g).first;
        edge_color[e] = 1;

        alps::graph::canonical_properties_type<graph_type>::type gp(canonical_properties(g));
        expect_partition(get<alps::graph::partition>(gp), {{0, 3, 5, 6}, {1, 2, 4, 7}});
    }

    //
    //  7---0---6---2       // c0 ---
    //  +           +       // c1 +++
    //  1---5---3---4
    //
    {
        graph_type g(8);
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 7, g).first;
        edge_color[e] = 1;
        e = add_edge(0, 6, g).first;
        edge_color[e] = 1;
        e = add_edge(2, 6, g).first;
        edge_color[e] = 1;
        e = add_edge(2, 4, g).first;
        edge_color[e] = 0;
        e = add_edge(3, 4, g).first;
        edge_color[e] = 1;
        e = add_edge(3, 5, g).first;
        edge_color[e] = 1;
        e = add_edge(1, 5, g).first;
        edge_color[e] = 1;
        e = add_edge(1, 7, g).first;
        edge_color[e] = 0;

        alps::graph::canonical_properties_type<graph_type>::type gp(canonical_properties(g));
        expect_partition(get<alps::graph::partition>(gp), {{0, 3, 5, 6}, {1, 2, 4, 7}});
    }
}
