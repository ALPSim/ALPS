/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013 by Andreas Hehn <hehn@phys.ethz.ch>                          *
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


#ifndef ALPS_GRAPH_UTILS_HPP
#define ALPS_GRAPH_UTILS_HPP

#include <alps/lattice/graphproperties.h>
#include <alps/lattice/propertymap.h>
#include <alps/graph/canonical_properties_traits.hpp>
#include <boost/static_assert.hpp>
#include <cassert>
#include <algorithm>
#include <vector>
#include <map>

namespace alps {
namespace graph {
namespace detail {
    template <typename Graph, typename PropertyTag>
    struct iteration_selector;

    template <typename Graph>
    struct iteration_selector<Graph,alps::vertex_type_t>
    {
        typedef typename boost::graph_traits<Graph>::vertex_iterator iterator;
        static std::pair<iterator,iterator> range(Graph const& g) { return vertices(g); }
    };

    template <typename Graph>
    struct iteration_selector<Graph,alps::edge_type_t>
    {
        typedef typename boost::graph_traits<Graph>::edge_iterator iterator;
        static std::pair<iterator,iterator> range(Graph const& g) { return edges(g); }
    };

    struct less_first
    {
        template <typename T>
        bool operator() (T const& a, T const& b)
        {
            return a.first < b.first;
        }
    };
    struct less_second
    {
        template <typename T>
        bool operator() (T const& a, T const& b)
        {
            return a.second < b.second;
        }
    };

    template <typename T2>
    struct equal_second
    {
        equal_second(T2 const& t)
        : v_(t)
        {
        }

        template <typename T1>
        bool operator() (std::pair<T1,T2> const& a) const
        {
            return a.second == v_;
        }

      private:
        T2 v_;
    };
}

/**
  * Extracts a sorted list of all vertex or edge colors in a graph.
  * \parm tag Selects if vertex or edge colors should be listed, it is either alps::vertex_type_t or alps::edge_type_t.
  * \parm g The graph to be inspected.
  * \return A vector containing allcolors that occur in the graph in ascending order.
  */
template <typename Graph, typename PropertyTag>
std::vector<typename boost::property_map<Graph,PropertyTag>::type::value_type> get_color_list(PropertyTag tag, Graph const& g)
{
    BOOST_STATIC_ASSERT(( boost::is_same<PropertyTag,alps::edge_type_t>::value || boost::is_same<PropertyTag,alps::vertex_type_t>::value ));
    BOOST_STATIC_ASSERT(( has_property<PropertyTag,Graph>::edge_property || has_property<PropertyTag,Graph>::vertex_property ));
    using std::sort;
    std::vector<typename boost::property_map<Graph,PropertyTag>::type::value_type> colors;
    typename boost::graph_traits<Graph>::edge_iterator e_it, e_end;

    typename detail::iteration_selector<Graph,PropertyTag>::iterator it,end;

    // We expect only a small number of different types
    for(boost::tie(it,end) = detail::iteration_selector<Graph,PropertyTag>::range(g); it != end; ++it)
    {
        typename boost::property_map<Graph,PropertyTag>::type::value_type ep = get(tag, g)[*it];
        if( find(colors.begin(),colors.end(),ep) == colors.end() )
            colors.push_back(ep);
    }
    sort(colors.begin(),colors.end());
    return colors;
}


/**
  * Checks if the color map has an entry for each edge color present in graph g.
  * \param g Graph this map will be applied to.
  * \param map A map that maps an edge type/color (unsigned int) of the graph to a new edge type/color (unsigned int).
  */
template <typename Graph>
bool color_map_covers_all_colors(Graph const& g, std::map<unsigned int, unsigned int> const& map)
{
    typedef typename boost::property_map<Graph, edge_type_t>::type::value_type edge_type;
    BOOST_STATIC_ASSERT(( boost::is_same< edge_type, unsigned int>::value ));
    std::vector<edge_type> colorlist = get_color_list(alps::edge_type_t(),g); 
    for(typename std::vector<edge_type>::const_iterator it=colorlist.begin(); it != colorlist.end(); ++it)
    {
        if(map.count(*it) == 0)
            return false;
    }
    return true;
}

/**
  * Remaps the edge types of graph g according to the specified map.
  * This function modifies the edge properties of the grap.
  * \param g Graph to modify.
  * \param map A map that maps an edge type (unsigned int) of the graph to a new edge type (unsigned int).
  */
template <typename Graph>
void remap_edge_types(Graph& g, std::map<unsigned int, unsigned int> const& map)
{
    typedef typename boost::property_map<Graph, edge_type_t>::type::value_type edge_type;
    BOOST_STATIC_ASSERT(( boost::is_same< edge_type, unsigned int>::value ));
    if(!color_map_covers_all_colors(g, map))
        throw std::runtime_error("Map does not cover all colors present in the graph.");
    typename boost::graph_traits<Graph>::edge_iterator it, end;
    for(boost::tie(it,end) = edges(g); it != end; ++it)
    {
        unsigned int type = get(alps::edge_type_t(),g,*it);
        unsigned int mappedtype = map.find(type)->second; // fine after color_map_covers_all_colors
        put(alps::edge_type_t(),g,*it,mappedtype);
    }
}

/**
  * Constructs all color maps which map to colors that are equivalent under the symmetry given by the color_partition
  * \parm The graph type, to deduce the color partition type
  * \parm color_partition the color partition defining the symmetries between the colors.
  * \return a vector of color mappings (which is a vector itself)
  */
template <typename Graph>
std::vector< std::vector<typename boost::property_map<Graph,alps::edge_type_t>::type::value_type> >
get_all_color_mappings_from_color_partition(Graph const&, typename color_partition<Graph>::type const& color_partition)
{
    using std::sort;
    using std::next_permutation;
    using std::max_element;
    typedef typename boost::property_map<Graph,alps::edge_type_t>::type::value_type color_type;
    unsigned int const end_c_id  = max_element(color_partition.begin(),color_partition.end(), detail::less_first())->first + 1;
    unsigned int const end_cg_id = max_element(color_partition.begin(),color_partition.end(), detail::less_second())->second + 1;

    std::vector< std::vector<color_type> > result( 1, std::vector<color_type>() );
    std::vector< std::vector<color_type> > colors_per_partition(end_cg_id, std::vector<color_type>());
    for(color_type c=0; c < end_c_id; ++c)
    {
        result[0].push_back(c);
        colors_per_partition[color_partition.find(c)->second].push_back(c);
    }

    typename std::vector< std::vector<color_type> >::iterator cg_it(colors_per_partition.begin()), cg_end(colors_per_partition.end());
    std::vector<color_type> permutations;
    for(; cg_it != cg_end; ++cg_it)
    {
        permutations.clear();
        permutations.insert(permutations.end(), cg_it->begin(), cg_it->end());
        sort(permutations.begin(), permutations.end());
        std::size_t const old_result_size = result.size();
        // We skip the first permutation, because that's the identity.
        while (next_permutation(permutations.begin(),permutations.end()))
        {
            result.reserve(result.size() + old_result_size);
            typename std::vector< std::vector<color_type> >::iterator rit(result.begin()), rend(result.begin() + old_result_size);
            for(; rit != rend; ++rit)
            {
                result.push_back(*rit);
                for(std::size_t i=0; i < permutations.size(); ++i)
                    result.back()[(*cg_it)[i]] = permutations[i];
            }
        }
    }
    return result;
}

} // end namespace graph
} // end namespace alps

#endif //ALPS_GRAPH_UTILS_HPP
