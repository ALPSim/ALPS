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
  
BOOST_AUTO_TEST_CASE(test_ctor_in_modular_accum)
{
    typedef alps::accumulator::accumulator<int
                                    , alps::accumulator::features<
                                          alps::accumulator::tag::fixed_size_binning
                                        , alps::accumulator::tag::max_num_binning
                                        , alps::accumulator::tag::autocorrelation>
                                    > accum;
    accum acc(alps::accumulator::bin_num = 10, alps::accumulator::bin_size = 10);
    
    acc << 1;
    acc << 2;
    acc << 3;
    acc << 4;
    acc << 5;
    //~ 
    accum acc2(acc);
        //~ 
    BOOST_REQUIRE(count(acc2) == acc.count());
    BOOST_REQUIRE( mean(acc2) == acc.mean());
    BOOST_REQUIRE( error(acc2) == acc.error());
    BOOST_REQUIRE( fixed_size_binning(acc2).bins() == acc.fixed_size_binning().bins());
    BOOST_REQUIRE( max_num_binning(acc2).bins() == acc.max_num_binning().bins());
    BOOST_REQUIRE( autocorrelation(acc2).bins() == acc.autocorrelation().bins());
    
}
