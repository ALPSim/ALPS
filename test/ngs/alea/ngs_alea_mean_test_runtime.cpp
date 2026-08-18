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

BOOST_AUTO_TEST_CASE(test_mean_in_modular_accum)
{
    alps::accumulator::accumulator<int, alps::accumulator::features<alps::accumulator::tag::mean> > acci;
    
    for(int i = 0; i < 101; ++i)
        acci << i;
        
    BOOST_REQUIRE( alps::accumulator::mean(acci) == 50);
    
    
    alps::accumulator::accumulator<double, alps::accumulator::features<alps::accumulator::tag::mean> > accd;
    
    for(double i = 0; i < 1.01; i += .01)
        accd << i;
        
    BOOST_REQUIRE( alps::accumulator::mean(accd) > .49999999999);
    BOOST_REQUIRE( alps::accumulator::mean(accd) < .50000000001);
}
