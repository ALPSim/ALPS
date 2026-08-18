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

#include <iostream>

#include "mean_archetype.hpp"

BOOST_AUTO_TEST_CASE(test_stream_with_mean_archetype)
{
    typedef alps::accumulator::accumulator<mean_archetype, alps::accumulator::features<alps::accumulator::tag::mean> > accum;
    accum acc;
    alps::accumulator::detail::accumulator_wrapper m(acc);
    
    for(int i = 0; i < 10; ++i)
        acc << mean_archetype();
            
    for(int i = 0; i < 10; ++i)
        m << mean_archetype();
        
    std::cout << acc << std::endl;
    std::cout << m << std::endl;    
}
