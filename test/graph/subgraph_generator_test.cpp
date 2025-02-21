/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                              Lukas Gamper <gamperl@gmail.com>                   *
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

//#define USE_COMPRESSED_EMBEDDING2

#include <alps/graph/subgraph_generator.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <vector>
#include <iostream>



enum { test_graph_size = 10 };

template <typename Graph>
void subgraph_generator_test(unsigned int order_ )
{
    std::ifstream in("../../lib/xml/lattices.xml");
    alps::Parameters parm;
    parm["LATTICE"] = "square lattice";
    parm["L"] = 2*order_+1;
    alps::graph_helper<> alps_lattice(in,parm);

    typedef alps::coordinate_graph_type lattice_graph_type;
    lattice_graph_type lattice_graph = alps_lattice.graph();


    typedef alps::graph::subgraph_generator<Graph,lattice_graph_type> graph_gen_type;
    std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> pin(1, 2*order_*order_);
    graph_gen_type graph_gen(lattice_graph, pin);

    typename graph_gen_type::iterator it,end;
    boost::tie(it,end) = graph_gen.generate_up_to_n_edges(order_);
    std::cout<<std::distance(it,end)<<std::endl;
}

int main()
{
    typedef boost::adjacency_list<boost::vecS, boost::vecS,boost::undirectedS> graph_type;
    subgraph_generator_test<graph_type>(test_graph_size);
    return 0;
}
