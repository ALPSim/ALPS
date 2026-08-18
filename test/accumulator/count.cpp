/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/config.hpp>

#ifndef ALPS_NGS_USE_NEW_ALEA
#error "This test only works with new alea library"
#endif

#define BOOST_TEST_MODULE alps::ngs::accumulator

#include <alps/ngs/accumulator.hpp>

#ifndef ALPS_LINK_BOOST_TEST
#	include <boost/test/included/unit_test.hpp>
#else
#	include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE(count_feature) {

	alps::accumulator::accumulator_set measurements;
	measurements << alps::accumulator::RealObservable("scalar")
				 << alps::accumulator::RealVectorObservable("vector");

	for (int i = 1; i < 1001; ++i) {
		measurements["scalar"] << i;
		BOOST_REQUIRE(count(measurements["scalar"]) == i);
		measurements["vector"] << std::vector<double>(10, i);
		BOOST_REQUIRE(count(measurements["vector"]) == i);
	}

	alps::accumulator::result_set results(measurements);
	BOOST_REQUIRE(count(results["scalar"]) == 1000);
	BOOST_REQUIRE(count(results["vector"]) == 1000);
}
