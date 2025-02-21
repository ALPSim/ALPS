/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
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

BOOST_AUTO_TEST_CASE(mean_feature) {

	alps::accumulator::accumulator_set measurements;
	measurements << alps::accumulator::RealObservable("obs1")
				 << alps::accumulator::RealObservable("obs2");

	for (int i = 1; i < 1000; ++i) {
		measurements["obs1"] << 1.;
		BOOST_REQUIRE(measurements["obs1"].mean<double>() == 1.);
		measurements["obs2"] << i;
		BOOST_REQUIRE(measurements["obs2"].mean<double>() == double(i + 1) / 2.);
	}

	alps::accumulator::result_set results(measurements);
	BOOST_REQUIRE(results["obs1"].mean<double>() == 1.);
	BOOST_REQUIRE(results["obs2"].mean<double>() == 500.);
}
