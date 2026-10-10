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

BOOST_AUTO_TEST_CASE(test_count_in_modular_accum)
{
    alps::accumulator::accumulator<int> acci;
    
    for (int i = 0; i < 100; ++i)
        acci << i;
        
    BOOST_REQUIRE( alps::accumulator::count(acci) == 100);
    
    
    alps::accumulator::accumulator<double> accd;
    
    for (double i = 0; i < 1; i += .01)
        accd << i;
            
    BOOST_REQUIRE( alps::accumulator::count(accd) == 100);
}
