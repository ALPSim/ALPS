/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013 by Andreas Hehn <hehn@phys.ethz.ch>                          *
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

template <typename Partition>
void dump_partition(std::ostream& os, Partition const& partition)
{
    os << "(";
    for (typename Partition::const_iterator it = partition.begin(); it != partition.end(); ++it) {
        os << "(";
        for (typename Partition::value_type::const_iterator jt = it->begin(); jt != it->end(); ++jt)
            os << (jt == it->begin() ? "" : " ") << *jt;
        os << ")";
    }
    os << ")";
}

void orbit_test1()
{
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

    dump_partition(std::cout, get<alps::graph::partition>(gp));
    std::cout << std::endl;
    dump_partition(std::cout, get<alps::graph::partition>(gp_with_sym));
    std::cout << std::endl;
}

void orbit_test2()
{
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

    dump_partition(std::cout, get<alps::graph::partition>(gp));
    std::cout << std::endl;
    dump_partition(std::cout, get<alps::graph::partition>(gp_with_sym));
    std::cout << std::endl;
}

void orbit_test3()
{
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

    dump_partition(std::cout, get<alps::graph::partition>(gp));
    std::cout << std::endl;
    dump_partition(std::cout, get<alps::graph::partition>(gp_with_sym));
    std::cout << std::endl;
}

void orbit_test4()
{
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
        dump_partition(std::cout, get<alps::graph::partition>(gp));
        std::cout << std::endl;
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
        dump_partition(std::cout, get<alps::graph::partition>(gp));
        std::cout << std::endl;
    }
}

void orbit_test5()
{
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
        dump_partition(std::cout, get<alps::graph::partition>(gp));
        std::cout << std::endl;
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
        dump_partition(std::cout, get<alps::graph::partition>(gp));
        std::cout << std::endl;
    }
}
int main() {
    orbit_test1();
    orbit_test2();
    orbit_test3();
    orbit_test4();
    orbit_test5();
    return 0;
}
