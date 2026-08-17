/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2016 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/graph/lattice_constant.hpp>

#include <boost/progress.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <iostream>

int main() {
    using boost::get;
    using alps::graph::canonical_properties;

    typedef unsigned int lc_type;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_type;

    alps::Parameters parm;
    unsigned int side_length = 40;

    std::ifstream in("../../lib/xml/lattices.xml");
    parm["LATTICE"] = "triangular lattice";
    parm["L"] = side_length;
    alps::graph_helper<> lattice(in,parm);

    std::vector<std::pair<graph_type,lc_type> > g;


    //
    //  0---1
    //  |   |
    //  2---3
    //
    g.push_back(std::make_pair(graph_type(), 3)); // verified on paper
    add_edge(0, 1, g.back().first);
    add_edge(1, 3, g.back().first);
    add_edge(3, 2, g.back().first);
    add_edge(2, 0, g.back().first);

    //
    //  1---0---2
    //
    g.push_back(std::make_pair(graph_type(),15)); // verified on paper
    add_edge(0, 1,g.back().first);
    add_edge(0, 2,g.back().first);

    //
    //  3---1---0---2
    //
    g.push_back(std::make_pair(graph_type(),69)); // verified on paper
    add_edge(0, 1,g.back().first);
    add_edge(0, 2,g.back().first);
    add_edge(1, 3,g.back().first);

    //
    //  1---0---2
    //      |
    //      3
    //
    g.push_back(std::make_pair(graph_type(),20)); // verified on paper
    add_edge(0, 1, g.back().first);
    add_edge(0, 2, g.back().first);
    add_edge(0, 3, g.back().first);


    //
    //     3
    //     |
    // 1---0---2
    //     |
    //     4
    //
    g.push_back(std::make_pair(graph_type(),15)); // verified on paper
    add_edge(0, 1,g.back().first);
    add_edge(0, 2,g.back().first);
    add_edge(0, 3,g.back().first);
    add_edge(0, 4,g.back().first);

    //
    //   2   5
    //   |   |
    //   0---1
    //   |   |
    //   3   4
    //
    g.push_back(std::make_pair(graph_type(),207)); // just saved output of hopefully correct impl
    add_edge(0, 1,g.back().first);
    add_edge(0, 2,g.back().first);
    add_edge(0, 3,g.back().first);
    add_edge(1, 4,g.back().first);
    add_edge(1, 5,g.back().first);

    //
    //           8
    //           |
    //           4
    //           |
    //   6---2---0---1---5
    //           |
    //           3
    //           |
    //           7
    //
    g.push_back(std::make_pair(graph_type(),2280)); // just saved output of hopefully correct impl
    add_edge(0, 1,g.back().first);
    add_edge(0, 2,g.back().first);
    add_edge(0, 3,g.back().first);
    add_edge(0, 4,g.back().first);
    add_edge(1, 5,g.back().first);
    add_edge(2, 6,g.back().first);
    add_edge(3, 7,g.back().first);
    add_edge(4, 8,g.back().first);

    //
    //         3
    //         |
    //   6--2--0--1--4--7
    //            |
    //            5
    //
    g.push_back(std::make_pair(graph_type(),12438)); // just saved output of hopefully correct impl
    add_edge(0, 1,g.back().first);
    add_edge(0, 2,g.back().first);
    add_edge(0, 3,g.back().first);
    add_edge(1, 4,g.back().first);
    add_edge(1, 5,g.back().first);
    add_edge(2, 6,g.back().first);
    add_edge(4, 7,g.back().first);

    int success = 0;
    {
        boost::progress_timer timer;
        for(std::vector<std::pair<graph_type,lc_type> >::iterator it= g.begin(); it != g.end(); ++it)
        {
            lc_type lc = alps::graph::lattice_constant(
                  it->first
                , lattice.graph()
                , lattice.lattice()
                , alps::cell(std::vector<int>(2,side_length/2),lattice.lattice()) //side_length * side_length / 2 + side_length / 2 - 1
            );
            if (lc != it->second) {
                std::cerr<<"ERROR: lattice constant does not match!"<<std::endl;
                std::cerr<<"Calculated: "<<lc<<"\tReference: "<<it->second<<std::endl<<std::endl;
                success = -1;
            } else
                std::cerr<<"SUCCESS: "<<it->second<<std::endl;
        }
    }
    return success;
}
