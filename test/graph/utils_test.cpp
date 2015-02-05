/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2015 by Andreas Hehn <hehn@phys.ethz.ch>                          *
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
#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include <alps/graph/utils.hpp>
#include <boost/container/flat_map.hpp>

typedef boost::property<alps::edge_type_t,unsigned int> edge_props;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, edge_props> graph_type;

bool get_all_color_mappings_from_color_partition_test()
{
    typedef std::vector< std::vector<alps::type_type> > color_mapping_type;

    alps::graph::color_partition<graph_type>::type color_symmetry;
    graph_type g;
    color_symmetry[0] = 0;
    color_symmetry[1] = 1;
    color_symmetry[2] = 0;
    color_symmetry[3] = 0;
    color_symmetry[4] = 1;
    color_symmetry[5] = 1;
    color_symmetry[6] = 2;
    color_symmetry[7] = 0;
    color_symmetry[8] = 0;
    
    color_mapping_type color_mappings(alps::graph::get_all_color_mappings_from_color_partition(g, color_symmetry));

    for(color_mapping_type::const_iterator it = color_mappings.begin(); it != color_mappings.end(); ++it)
    {
        for(std::size_t c = 0; c < it->size(); ++c)
            std::cout << c << "->" << (*it)[c] << ", ";
        std::cout << std::endl;
    }
    std::cout << "Num. of color_mappings: " << color_mappings.size() << std::endl;
    return color_mappings.size() == 720; //= 5! * 3! * 1!
}

int main()
{
    bool ok = get_all_color_mappings_from_color_partition_test();
    if(!ok)
        std::cout << "get_all_color_mappings_from_color_partition_test FAILED" << std::endl;
    return ok ? 0 : 1;
}
