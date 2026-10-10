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

template<typename T>
alps::accumulator::autocorrelation_proxy_type<T> fct()
{
    return alps::accumulator::autocorrelation_proxy_type<T>();
}

BOOST_AUTO_TEST_CASE(test_wrapper_for_modular_accum)
{
    fct<int>();
    //~ typedef alps::accumulator::accumulator<int, alps::accumulator::features<alps::accumulator::tag::mean> > accum;
    //~ accum acci;
    //~ 
    //~ alps::accumulator::detail::accumulator_wrapper m(acci);
    //~ 
    //~ for(int i = 0; i < 101; ++i)
    //~ {
        //~ m << i;
    //~ }
        
    //~ BOOST_REQUIRE( m.get<int>().mean() == 50);
    //~ BOOST_REQUIRE( alps::accumulator::mean(alps::accumulator::extract<accum>(m)) == 50);
    //~ BOOST_REQUIRE( alps::accumulator::mean(m.extract<accum>()) == 50);
    
}
