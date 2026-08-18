/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2015 by Andreas Hehn <hehn@phys.ethz.ch>                          *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_GRAPH_DETAIL_HELPER_FUNCTIONS_HPP
#define ALPS_GRAPH_DETAIL_HELPER_FUNCTIONS_HPP

#include <alps/graph/canonical_properties_traits.hpp>
#include <vector>
#include <boost/type_traits/is_unsigned.hpp>
#include <boost/static_assert.hpp>

namespace alps {
namespace graph {
namespace detail {

template <typename Graph>
std::vector<std::size_t> build_vertex_to_partition_index_map(typename alps::graph::partition_type<Graph>::type const& partition, Graph const& G)
{
    BOOST_STATIC_ASSERT(( boost::is_unsigned<typename boost::graph_traits<Graph>::vertex_descriptor>::value ));
    std::vector<std::size_t> partition_of(num_vertices(G));
    for (typename alps::graph::partition_type<Graph>::type::const_iterator it = partition.begin(); it != partition.end(); ++it)
    {
        for (typename alps::graph::partition_type<Graph>::type::value_type::const_iterator jt = it->begin(); jt != it->end(); ++jt)
        {
            assert( *jt < num_vertices(G) && "vertex descriptors need to be unsigned numbers in range [0,num_vertices(G)[");
            partition_of[*jt] = it - partition.begin();
        }
    }
    return partition_of;
}

} // end namespace detail
} // end namespace graph
} // end namespace alps


#endif // ALPS_GRAPH_DETAIL_HELPER_FUNCTIONS_HPP
