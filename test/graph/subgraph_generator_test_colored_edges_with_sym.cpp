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
