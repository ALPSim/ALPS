/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2012 by Andreas Hehn <hehn@phys.ethz.ch>                          *
 *                       Lukas Gamper <gamperl@gmail.com>                          *
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

#include <alps/lattice.h>
#include <alps/graph/lattice_constant.hpp>

#include <boost/graph/adjacency_list.hpp>

#include <iostream>

int main() {
    using boost::get;
    using boost::put;
    using alps::graph::canonical_properties;
    using alps::graph::canonical_properties_type;

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_type;

    typedef alps::graph_helper<>::vertex_iterator vertex_iterator;
    typedef alps::graph_helper<>::edge_iterator edge_iterator;
    
    alps::Parameters parm;
    unsigned int side_length = 40;
    
    std::ifstream in("../../lib/xml/lattices.xml");
    parm["LATTICE"] = "square lattice";
    parm["L"] = side_length;
    alps::graph_helper<>::vertex_descriptor center_vertex = side_length * side_length / 2 + side_length / 2;

    alps::graph_helper<> lattice(in,parm);
    
    graph_type lattice_graph(num_vertices(lattice.graph()));
    boost::graph_traits<alps::graph_helper<>::graph_type>::edge_iterator it, et;
    for(boost::tie(it, et) = edges(lattice.graph()); it != et; ++it)
    {
        add_edge(source(*it, lattice.graph()), target(*it, lattice.graph()), lattice_graph);
    }

    typedef unsigned int contrib_type;
    typedef boost::tuple<canonical_properties_type<graph_type>::type, unsigned int, alps::numeric::vector<contrib_type> > input_type;
    typedef alps::numeric::vector<contrib_type> output_type;
    std::vector<boost::tuple<graph_type, input_type, output_type> > test_graphs;


    // An initial output vector filled with 'random' data
    output_type init(num_vertices(lattice_graph));
    for(std::size_t i=0; i < init.size();++i)
        init[i] = i;

    //
    // 1---0
    //
    {
        graph_type g;
        add_edge(0,1,g);

        // Orbit partition: (0 1) -> [0] (1)
        unsigned int breaking_vertex = 0;
        alps::numeric::vector<contrib_type> part_contrib(2);
        part_contrib[0] = 2; // c[0 -> 0]
        part_contrib[1] = 3; // c[0 -> 1]
        input_type in(canonical_properties(g,breaking_vertex), breaking_vertex, part_contrib);

        output_type out(init);
        out[820] += 2*2*2;  // (0,0)
        out[821] += 1*3;    // (1,0)
        out[819] += 1*3;    // (-1,0)
        out[860] += 1*3;    // (0,1)
        out[780] += 1*3;    // (0,-1)

        test_graphs.push_back(boost::make_tuple(g,in,out));
    }

    //
    // 1---0---2
    //     |
    //     3
    //

    {
        graph_type g;
        add_edge(0,1,g);
        add_edge(0,2,g);
        add_edge(0,3,g);

        {
            // Orbit partition: (0)(1 2 3) -> [0] (1 2 3)
            unsigned int breaking_vertex = 0;
            alps::numeric::vector<contrib_type> part_contrib(2);
            part_contrib[0] = 2;
            part_contrib[1] = 3;
            input_type in(canonical_properties(g,breaking_vertex), breaking_vertex, part_contrib);

            output_type out(init);
            out[820] += 4*2;  // (0,0)
            out[821] += 3*3;  // (1,0)
            out[819] += 3*3;  // (-1,0)
            out[860] += 3*3;  // (0,1)
            out[780] += 3*3;  // (0,-1)

            test_graphs.push_back(boost::make_tuple(g,in,out));
        }
        {
            // Orbit partition: (0)(1 2 3) -> (0) [1] (2 3)
            unsigned int breaking_vertex = 1;
            alps::numeric::vector<contrib_type> part_contrib(3);
            part_contrib[0] = 5;  // c[1->0]
            part_contrib[1] = 7;  // c[1->1]
            part_contrib[2] = 11; // c[1->2]
            input_type in(canonical_properties(g,breaking_vertex), breaking_vertex, part_contrib);

            output_type out(init);
            out[820]  += 12 * part_contrib[1];  // ( 0, 0)
            out[821]  += 3  * part_contrib[0];  // ( 1, 0)
            out[860]  += 3  * part_contrib[0];  // ( 0, 1)
            out[819]  += 3  * part_contrib[0];  // (-1, 0)
            out[780]  += 3  * part_contrib[0];  // ( 0,-1)

            out[861] += 4*11; // ( 1, 1)
            out[859] += 4*11; // (-1, 1)
            out[781] += 4*11; // ( 1,-1)
            out[779] += 4*11; // (-1,-1)

            out[822] += 2*11; // ( 2, 0)
            out[818] += 2*11; // (-2, 0)
            out[900] += 2*11; // ( 0, 2)
            out[740] += 2*11; // ( 0,-2)

            test_graphs.push_back(boost::make_tuple(g,in,out));
        }
    }

    int success = 0;
    for(std::vector<boost::tuple<graph_type, input_type, output_type> >::iterator it = test_graphs.begin(); it != test_graphs.end(); ++it) {

        alps::numeric::matrix<unsigned int> embedding_factors(num_vertices(lattice_graph), get<2>(get<1>(*it)).size());

        output_type output(init);
        alps::graph::lattice_constant(
              embedding_factors
            , get<0>(*it)
            , lattice_graph
            , lattice.lattice()
            , center_vertex
            , get<1>(get<1>(*it))
            , get<alps::graph::partition>(get<0>(get<1>(*it)))
        );
        output += embedding_factors * get<2>(get<1>(*it));
        output_type ref = get<2>(*it);
        if ( output != ref )
        {
            std::cerr<<"ERROR: lattice constant does not match!"<<std::endl;
            std::cerr<<"Graph:"<<std::distance(test_graphs.begin(),it)<<":"<<std::endl;
            std::cerr<<"Calculated <-> Reference"<<std::endl;
            for(std::size_t i=0; i != output.size(); ++i)
            {
                if(output[i] != ref[i])
                    std::cerr<<"vertex: "<<i<<":\t"<<output[i]-init[i]<<"\t<->\t"<<ref[i]-init[i]<<std::endl;
            }
            std::cerr<<std::endl;
            success = -1;
        }
    }
    return success;
}
