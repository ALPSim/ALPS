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
#include <gtest/gtest.h>
#include "graph_assertions.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include <set>
#include <alps/graph/utils.hpp>
#include <boost/container/flat_map.hpp>

typedef boost::property<alps::edge_type_t,unsigned int> edge_props;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, edge_props> graph_type;

TEST(GraphUtilities, EnumeratesAllColorPermutations)
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

    EXPECT_EQ(color_mappings.size(), 720u); // 5! * 3! * 1!
    const std::set<std::vector<alps::type_type>> unique(color_mappings.begin(), color_mappings.end());
    EXPECT_EQ(unique.size(), color_mappings.size());
    for (const auto& mapping : color_mappings) {
        ASSERT_EQ(mapping.size(), color_symmetry.size());
        expect_vertex_permutation(mapping, color_symmetry.size());
        for (std::size_t color = 0; color < mapping.size(); ++color)
            EXPECT_EQ(color_symmetry[color], color_symmetry[mapping[color]]);
    }
}
