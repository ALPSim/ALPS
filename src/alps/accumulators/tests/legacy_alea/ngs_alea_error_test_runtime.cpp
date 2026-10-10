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
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_CASE(test_error_in_modular_accum)
{
    alps::accumulator::accumulator<int, alps::accumulator::features<alps::accumulator::tag::error> > acci;
    
    acci << 2;
    acci << 6;
    
    BOOST_REQUIRE( error(acci) == 2.);
    
    
    alps::accumulator::accumulator<double, alps::accumulator::features<alps::accumulator::tag::error> > accd;
    
    
    accd << .2;
    accd << .6;
    
    BOOST_REQUIRE_CLOSE(alps::accumulator::error(accd), 0.2, 0.01);
}
