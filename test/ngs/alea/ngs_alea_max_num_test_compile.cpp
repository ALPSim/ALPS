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

#include "error_archetype.hpp"
 
BOOST_AUTO_TEST_CASE(test_max_num_bin_with_error_archetype)
{
    alps::accumulator::accumulator<error_archetype, alps::accumulator::features<alps::accumulator::tag::max_num_binning> > acci; //Default 128
    
    acci << 2;
    acci << 6;
    
    BOOST_REQUIRE( alps::accumulator::max_num_binning(acci).bin_number() == 128);
    
    alps::accumulator::accumulator<error_archetype, alps::accumulator::features<alps::accumulator::tag::max_num_binning> > accd(alps::accumulator::bin_num = 10);
    
    accd << .2;
    accd << .6;
        
    BOOST_REQUIRE( alps::accumulator::max_num_binning(accd).bin_number() == 10);
    
    alps::accumulator::max_num_binning(accd);
    accd.max_num_binning();
}
