/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2012 - 2015 by Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                              Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/graph/lattice_constant.hpp>
#include <alps/lattice/lattice.h>
#include <boost/graph/adjacency_list.hpp>
#include <iostream>

template <typename AlpsLattice>
typename alps::lattice_traits<AlpsLattice>::cell_descriptor get_middle_cell(AlpsLattice const& lattice)
{
    typename alps::lattice_traits<AlpsLattice>::offset_type offset = lattice.extent();
    for(unsigned int d=0; d < dimension(lattice); ++d)
        offset[d] /= 2;
    return cell(offset,lattice);
}

int main() {
    using boost::get;
    using boost::put;
    using alps::graph::canonical_properties;

    typedef unsigned int lc_type;
    typedef boost::property<alps::edge_type_t,alps::type_type> edge_props;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,boost::no_property,edge_props> graph_type;
    typedef alps::graph_helper<>::edge_iterator edge_iterator;
    typedef boost::graph_traits<graph_type>::edge_descriptor  edge_descriptor;
    typedef boost::property_map<graph_type,alps::edge_type_t>::type edge_color_map_type;

    alps::Parameters parm;
    unsigned int side_length = 40;

    std::ifstream in("../../lib/xml/lattices.xml");
    parm["LATTICE"] = "anisotropic square lattice";
    parm["L"] = side_length;
    alps::graph_helper<> lattice(in,parm);

    std::vector<std::pair<graph_type,lc_type> > tests;

    // edge color 0 ...
    // edge color 1 ___
    
    //
    //  0...1
    //  |   |
    //  2...3
    //
    {
        graph_type g;
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        put(alps::edge_type_t(), g, e, 0);
        e = add_edge(1, 3, g).first;
        put(alps::edge_type_t(), g, e, 1);
        e = add_edge(3, 2, g).first;
        put(alps::edge_type_t(), g, e, 0);
        e = add_edge(2, 0, g).first;
        put(alps::edge_type_t(), g, e, 1);
        tests.push_back(std::make_pair(g, 1));
    }
    
    //
    //  0...1
    //  .   |
    //  2___3
    //
    {
        graph_type g;
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        put(alps::edge_type_t(), g, e, 0);
        e = add_edge(1, 3, g).first;
        put(alps::edge_type_t(), g, e, 1);
        e = add_edge(3, 2, g).first;
        put(alps::edge_type_t(), g, e, 1);
        e = add_edge(2, 0, g).first;
        put(alps::edge_type_t(), g, e, 0);
        tests.push_back(std::make_pair(g, 0));
    }



    //
    //  1___0___2
    //
    {
        graph_type g;
        edge_descriptor e;
        e = add_edge(0, 1,g).first;
        put(alps::edge_type_t(), g, e, 0);
        e = add_edge(0, 2,g).first;
        put(alps::edge_type_t(), g, e, 0);
        tests.push_back(std::make_pair(g,1));
    }
    
    //
    //  1...0___2
    //
    {
        graph_type g;
        edge_descriptor e;
        e = add_edge(0, 1,g).first;
        put(alps::edge_type_t(), g, e, 0);
        e = add_edge(0, 2,g).first;
        put(alps::edge_type_t(), g, e, 1);
        tests.push_back(std::make_pair(g,4));
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
        tests.push_back(std::make_pair(g,1));
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
        tests.push_back(std::make_pair(g,1));
    }

    int success = 0;

    for(std::vector<std::pair<graph_type, lc_type> >::iterator it = tests.begin(); it != tests.end(); ++it)
    {
        lc_type lc = alps::graph::lattice_constant(
              it->first
            , lattice.graph()
            , lattice.lattice()
            , get_middle_cell(lattice.lattice())
        );
        if (lc != it->second)
        {
            std::cerr<<"ERROR: lattice constant does not match!"<<std::endl;
            std::cerr<<"Graph:"<<std::distance(tests.begin(),it)<<" Calculated: "<<lc<<"\tReference: "<<it->second<<std::endl<<std::endl;
            success = -1;
        }
    }
    return success;
}
