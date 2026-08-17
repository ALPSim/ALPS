/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Mario Koenz <mkoenz@ethz.ch>                       *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define BOOST_TEST_MODULE alps::ngs::accumulator

#include <alps/ngs.hpp>

#ifndef ALPS_LINK_BOOST_TEST
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

#include "hist_archetype.hpp"

#include <utility>


BOOST_AUTO_TEST_CASE(test_hist_with_value_archetype)
{
    alps::accumulator::histogram_old<hist_archetype> acc(1, 6, 6);
    
    for(int i = 0; i < 10; ++i)
        acc << hist_archetype();
    acc[1] += hist_archetype();
    acc[1] ++;
    ++acc[1];
    acc << std::pair<int, hist_archetype>(1,hist_archetype());
    std::cout << acc[1] << std::endl;
    std::cout << acc << std::endl;
    
    acc.count();    
    acc.mean();
}
