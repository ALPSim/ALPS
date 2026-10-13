/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <gtest/gtest.h>
#include "graph_assertions.hpp"
#include <alps/graph/canonical_properties.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <iostream>


using boost::get;
using alps::graph::canonical_properties;


TEST(GraphCanonicalLabel, colored_edges_test) {
    typedef boost::property<alps::edge_type_t,unsigned int> edge_props;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, edge_props> graph_type;
    typedef boost::graph_traits<graph_type>::edge_descriptor edge_descriptor;
    typedef boost::property_map<graph_type,alps::edge_type_t>::type edge_color_map_type;

    graph_type g(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(1, 2, g).first;
        edge_color[e] = 0;
        e = add_edge(2, 3, g).first;
        edge_color[e] = 0;
        e = add_edge(0, 3, g).first;
        edge_color[e] = 0;
    }

    graph_type h(g);
    {
        edge_descriptor e = edge(2,3,h).first;
        get(alps::edge_type_t(),h)[e] = 1;
    }

    graph_type i(g);
    {
        edge_descriptor e = edge(2,3,i).first;
        get(alps::edge_type_t(),i)[e] = 2;
    }

    graph_type j(h);
    {
        edge_descriptor e = edge(2,3,j).first;
        get(alps::edge_type_t(),j)[e] = 1;
    }

    graph_type k(g);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),k);
        edge_descriptor e;
        e = edge( 0, 1, k).first;
        edge_color[e] = 1;
        e = edge( 1, 2, k).first;
        edge_color[e] = 1;
        e = edge( 2, 3, k).first;
        edge_color[e] = 1;
        e = edge( 0, 3, k).first;
        edge_color[e] = 1;
    }

    alps::graph::graph_label<graph_type>::type label_g(get<1>(canonical_properties(g)));
    alps::graph::graph_label<graph_type>::type label_h(get<1>(canonical_properties(h)));
    alps::graph::graph_label<graph_type>::type label_i(get<1>(canonical_properties(i)));
    alps::graph::graph_label<graph_type>::type label_j(get<1>(canonical_properties(j)));
    alps::graph::graph_label<graph_type>::type label_k(get<1>(canonical_properties(k)));

    expect_label_encoding(label_g, "(0001101100 1111 (0))");
    expect_label_encoding(label_h, "(0001101100 00011110 (0 1))");
    expect_label_encoding(label_i, "(0001101100 00011110 (0 2))");
    expect_label_encoding(label_j, "(0001101100 00011110 (0 1))");
    expect_label_encoding(label_k, "(0001101100 1111 (1))");

    EXPECT_FALSE((label_g == label_h));
    EXPECT_FALSE((label_h == label_i));
    EXPECT_TRUE((label_h == label_j));
    EXPECT_FALSE((label_g == label_k));


}

TEST(GraphCanonicalLabel, colored_edges_test2) {
    typedef boost::property<alps::edge_type_t,unsigned int> edge_props;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, edge_props> graph_type;
    typedef boost::graph_traits<graph_type>::edge_descriptor edge_descriptor;
    typedef boost::property_map<graph_type,alps::edge_type_t>::type edge_color_map_type;

    graph_type g;
    edge_descriptor e;
    e = add_edge(0,1,g).first;
    boost::put(alps::edge_type_t(),g,e,0);
    e = add_edge(0,2,g).first;
    boost::put(alps::edge_type_t(),g,e,0);

    graph_type h;
    e = add_edge(0,1,h).first;
    boost::put(alps::edge_type_t(),h,e,1);
    e = add_edge(0,2,h).first;
    boost::put(alps::edge_type_t(),h,e,1);

    graph_type i;
    e = add_edge(0,1,i).first;
    boost::put(alps::edge_type_t(),i,e,3);
    e = add_edge(0,2,i).first;
    boost::put(alps::edge_type_t(),i,e,3);

    graph_type j;
    e = add_edge(0,1,j).first;
    boost::put(alps::edge_type_t(),j,e,1);
    e = add_edge(0,2,j).first;
    boost::put(alps::edge_type_t(),j,e,3);

    graph_type k;
    e = add_edge(0,1,k).first;
    boost::put(alps::edge_type_t(),k,e,0);
    e = add_edge(0,2,k).first;
    boost::put(alps::edge_type_t(),k,e,1);


    alps::graph::graph_label<graph_type>::type label_g(get<1>(canonical_properties(g)));
    alps::graph::graph_label<graph_type>::type label_h(get<1>(canonical_properties(h)));
    alps::graph::graph_label<graph_type>::type label_i(get<1>(canonical_properties(i)));
    alps::graph::graph_label<graph_type>::type label_j(get<1>(canonical_properties(j)));
    alps::graph::graph_label<graph_type>::type label_k(get<1>(canonical_properties(k)));

    expect_label_encoding(label_g, "(010100 11 (0))");
    expect_label_encoding(label_h, "(010100 11 (1))");
    expect_label_encoding(label_i, "(010100 11 (3))");
    expect_label_encoding(label_j, "(010100 0110 (1 3))");
    expect_label_encoding(label_k, "(010100 0110 (0 1))");

    EXPECT_FALSE((label_g == label_h));
    EXPECT_FALSE((label_h == label_i));
    EXPECT_FALSE((label_h == label_j));
    EXPECT_FALSE((label_g == label_k));
}

TEST(GraphCanonicalLabel, simple_test) {
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_type;
    graph_type g;
    add_edge(0, 1,g);
    add_edge(0, 2,g);
    add_edge(0, 3,g);
    add_edge(2, 1,g);

    alps::graph::graph_label<graph_type>::type label(get<1>(canonical_properties(g)));

    expect_label_encoding(label, "(0101101000)");

}
