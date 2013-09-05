/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define BOOST_TEST_MODULE alps::ngs::accumulator

#include <alps/hdf5/array.hpp>
#include <alps/hdf5/multi_array.hpp>
#include <alps/ngs/alea/next/accumulator.hpp>

#ifndef ALPS_LINK_BOOST_TEST
#	include <boost/test/included/unit_test.hpp>
#else
#	include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE(ngs_alea_next) {
	using namespace alps::accumulator;

	accumulator_set::register_serializable_type<
		impl::Accumulator<double, count_tag, impl::AccumulatorBase<double> >
	>();

	accumulator_set::register_serializable_type<
		impl::Accumulator<double, mean_tag, impl::Accumulator<double, count_tag, impl::AccumulatorBase<double> > >
	>();

	accumulator_set::register_serializable_type<
		impl::Accumulator<double, error_tag, impl::Accumulator<double, mean_tag, impl::Accumulator<double, count_tag, impl::AccumulatorBase<double> > > >
	>();

	accumulator_set::register_serializable_type<
		impl::Accumulator<double, weight_holder_tag<
			impl::Accumulator<double, mean_tag, impl::Accumulator<double, count_tag, impl::AccumulatorBase<double> > >
		>, impl::Accumulator<double, mean_tag, impl::Accumulator<double, count_tag, impl::AccumulatorBase<double> > > >
	>();

	accumulator_set::register_serializable_type<
		impl::Accumulator<double, max_num_binning_tag, impl::Accumulator<
			double, error_tag, impl::Accumulator<double, mean_tag, impl::Accumulator<double, count_tag, impl::AccumulatorBase<double> > >
		> >
	>();

	accumulator_set::register_serializable_type<
		impl::Accumulator<std::vector<double>, mean_tag, impl::Accumulator<std::vector<double>, count_tag, impl::AccumulatorBase<std::vector<double> > > >
	>();

	accumulator_set::register_serializable_type<
		impl::Accumulator<std::vector<double>, max_num_binning_tag, impl::Accumulator<
			std::vector<double>, error_tag, impl::Accumulator<std::vector<double>, mean_tag, impl::Accumulator<std::vector<double>, count_tag, impl::AccumulatorBase<std::vector<double> > > >
		> >
	>();

	accumulator_set::register_serializable_type<
		impl::Accumulator<alps::multi_array<double, 3>, mean_tag, impl::Accumulator<alps::multi_array<double, 3>, count_tag, impl::AccumulatorBase<alps::multi_array<double, 3> > > >
	>();

	accumulator_set::register_serializable_type<
		impl::Accumulator<alps::multi_array<double, 3>, max_num_binning_tag, impl::Accumulator<
			alps::multi_array<double, 3>, error_tag, impl::Accumulator<alps::multi_array<double, 3>, mean_tag, impl::Accumulator<alps::multi_array<double, 3>, count_tag, impl::AccumulatorBase<alps::multi_array<double, 3> > > >
		> >
	>();

	accumulator_set accumulators;
	accumulators.insert("mean", boost::shared_ptr<accumulator_wrapper>(
		new accumulator_wrapper(impl::Accumulator<double, mean_tag, impl::Accumulator<double, count_tag, impl::AccumulatorBase<double> > >())
	));

	accumulators["mean"] << 1.;
	accumulators["mean"](2.);
	accumulators["mean"](4);
	accumulators["mean"](8l);
	accumulators["mean"](16.f);

    BOOST_REQUIRE(count(accumulators["mean"]) == 5);
    BOOST_REQUIRE(mean(accumulators["mean"].get<double>()) == 6.2);

    accumulators << SimpleRealObservable("error");

	accumulators["error"] << 6.;
	accumulators["error"](1.);
	accumulators["error"](1);
	accumulators["error"](1l);
	accumulators["error"](1.f);

    BOOST_REQUIRE(count(accumulators["error"]) == 5);
    BOOST_REQUIRE(mean(accumulators["error"].get<double>()) == 2);
    BOOST_REQUIRE(error(accumulators["error"].get<double>()) == 1);

	accumulators.insert("count", boost::shared_ptr<accumulator_wrapper>(new accumulator_wrapper(impl::Accumulator<double, count_tag, impl::AccumulatorBase<double> >())));
	for (int i = 0; i < 10; ++i)
		accumulators["count"] << 1.;

    BOOST_REQUIRE(count(accumulators["count"]) == 10);

	accumulators.insert("weighted", boost::shared_ptr<accumulator_wrapper>(
		new accumulator_wrapper(impl::Accumulator<double, weight_holder_tag<
			impl::Accumulator<double, mean_tag, impl::Accumulator<double, count_tag, impl::AccumulatorBase<double> > >
		>, impl::Accumulator<double, mean_tag, impl::Accumulator<double, count_tag, impl::AccumulatorBase<double> > > >())
	));

	accumulators["weighted"](1., 1.);
	accumulators["weighted"](2., 1.);
	accumulators["weighted"](4., 1.);
	accumulators["weighted"](8., 1.);
	accumulators["weighted"](16., 1.);

    BOOST_REQUIRE(count(accumulators["weighted"]) == 5);
    BOOST_REQUIRE(mean(accumulators["weighted"].get<double>()) == 6.2);

    BOOST_REQUIRE(count(*weight(accumulators["weighted"].get<double>())) == 5);
    BOOST_REQUIRE(mean(weight(accumulators["weighted"].get<double>())->get<double>()) == 1);

	accumulators.insert("vector", boost::shared_ptr<accumulator_wrapper>(
		new accumulator_wrapper(impl::Accumulator<std::vector<double>, mean_tag, impl::Accumulator<std::vector<double>, count_tag, impl::AccumulatorBase<std::vector<double> > > >())
	));

	accumulators["vector"](std::vector<double>(3, 1.));
	accumulators["vector"](std::vector<double>(3, 2.));
	accumulators["vector"](std::vector<double>(3, 4.));
	accumulators["vector"](std::vector<double>(3, 8.));
	accumulators["vector"](std::vector<double>(3, 16.));

    BOOST_REQUIRE(count(accumulators["vector"]) == 5);
    std::vector<double> vector_mean(3, 6.2);
    BOOST_REQUIRE(std::equal(vector_mean.begin(), vector_mean.end(), mean(accumulators["vector"].get<std::vector<double> >()).begin()));

	accumulators.insert("int", boost::shared_ptr<accumulator_wrapper>(
		new accumulator_wrapper(impl::Accumulator<int, mean_tag, impl::Accumulator<int, count_tag, impl::AccumulatorBase<int > > >())
	));
	accumulators["int"](1);

    BOOST_REQUIRE(count(accumulators["int"]) == 1);
    BOOST_REQUIRE(mean(accumulators["int"].get<int>()) == 1);

	accumulators.insert("vecint", boost::shared_ptr<accumulator_wrapper>(
		new accumulator_wrapper(impl::Accumulator<std::vector<int>, mean_tag, impl::Accumulator<std::vector<int>, count_tag, impl::AccumulatorBase<std::vector<int> > > >())
	));
	accumulators["vecint"](std::vector<int>(3, 1));

    BOOST_REQUIRE(count(accumulators["vecint"]) == 1);
    std::vector<int> vecint_mean(3, 1);
    BOOST_REQUIRE(std::equal(vecint_mean.begin(), vecint_mean.end(), mean(accumulators["vecint"].get<std::vector<int> >()).begin()));

	accumulators.insert("array", boost::shared_ptr<accumulator_wrapper>(
		new accumulator_wrapper(impl::Accumulator<boost::array<double, 3>, mean_tag, impl::Accumulator<boost::array<double, 3>, count_tag, impl::AccumulatorBase<boost::array<double, 3> > > >())
	));
	boost::array<double, 3> array_val = { {1., 2., 3.} };
	accumulators["array"](array_val);

    BOOST_REQUIRE(count(accumulators["array"]) == 1);
    BOOST_REQUIRE(std::equal(array_val.begin(), array_val.end(), mean(accumulators["array"].get<boost::array<double, 3> >()).begin()));

	accumulators.insert("multi_array", boost::shared_ptr<accumulator_wrapper>(
		new accumulator_wrapper(impl::Accumulator<alps::multi_array<double, 3>, mean_tag, impl::Accumulator<alps::multi_array<double, 3>, count_tag, impl::AccumulatorBase<alps::multi_array<double, 3> > > >())
	));
	alps::multi_array<double, 3> multi_array_val(boost::extents[2][2][2]);
	std::fill(multi_array_val.origin(), multi_array_val.origin() + 8, 1.);
	accumulators["multi_array"](multi_array_val);

    BOOST_REQUIRE(count(accumulators["multi_array"]) == 1);
    BOOST_REQUIRE(std::equal(multi_array_val.begin(), multi_array_val.end(), mean(accumulators["multi_array"].get<alps::multi_array<double, 3> >()).begin()));

	accumulators.insert("doublemaxbin", boost::shared_ptr<accumulator_wrapper>(
		new accumulator_wrapper(impl::Accumulator<double, max_num_binning_tag, impl::Accumulator<
			double, error_tag, impl::Accumulator<double, mean_tag, impl::Accumulator<double, count_tag, impl::AccumulatorBase<double > > >
		> >())
	));

	accumulators["doublemaxbin"] << 6.;
	accumulators["doublemaxbin"](1.);
	accumulators["doublemaxbin"](1);
	accumulators["doublemaxbin"](1l);
	accumulators["doublemaxbin"](1.f);

    BOOST_REQUIRE(count(accumulators["doublemaxbin"]) == 5);
    BOOST_REQUIRE(mean(accumulators["doublemaxbin"].get<double>()) == 2);
    BOOST_REQUIRE(error(accumulators["doublemaxbin"].get<double>()) == 1);
    // TODO: check binning ...

	accumulators.insert("vectormaxbin", boost::shared_ptr<accumulator_wrapper>(
		new accumulator_wrapper(impl::Accumulator<std::vector<double>, max_num_binning_tag, impl::Accumulator<
			std::vector<double>, error_tag, impl::Accumulator<std::vector<double>, mean_tag, impl::Accumulator<std::vector<double>, count_tag, impl::AccumulatorBase<std::vector<double> > > >
		> >())
	));

	accumulators["vectormaxbin"](std::vector<double>(3, 1.));
	accumulators["vectormaxbin"](std::vector<double>(3, 2.));
	accumulators["vectormaxbin"](std::vector<double>(3, 4.));
	accumulators["vectormaxbin"](std::vector<double>(3, 8.));
	accumulators["vectormaxbin"](std::vector<double>(3, 16.));

    BOOST_REQUIRE(count(accumulators["vectormaxbin"]) == 5);
    std::vector<double> vectormaxbin_mean(3, 6.2);
    BOOST_REQUIRE(std::equal(vectormaxbin_mean.begin(), vectormaxbin_mean.end(), mean(accumulators["vectormaxbin"].get<std::vector<double> >()).begin()));
    // TODO: check error ...
    // TODO: check binning ...

	accumulators.insert("multiarraymaxbin", boost::shared_ptr<accumulator_wrapper>(
		new accumulator_wrapper(impl::Accumulator<alps::multi_array<double, 3>, max_num_binning_tag, impl::Accumulator<
			alps::multi_array<double, 3>, error_tag, impl::Accumulator<alps::multi_array<double, 3>, mean_tag, impl::Accumulator<alps::multi_array<double, 3>, count_tag, impl::AccumulatorBase<alps::multi_array<double, 3> > > >
		> >())
	));

	alps::multi_array<double, 3> multi_array_max_bin_val(boost::extents[2][2][2]);
	std::fill(multi_array_max_bin_val.origin(), multi_array_max_bin_val.origin() + 8, 1.);
	accumulators["multiarraymaxbin"](multi_array_max_bin_val);
	accumulators["multiarraymaxbin"](multi_array_max_bin_val);
	accumulators["multiarraymaxbin"](multi_array_max_bin_val);
	accumulators["multiarraymaxbin"](multi_array_max_bin_val);
	accumulators["multiarraymaxbin"](multi_array_max_bin_val);

    BOOST_REQUIRE(count(accumulators["multiarraymaxbin"]) == 5);
    // TODO: check mean ...
    // TODO: check error ...
    // TODO: check binning ...

	std::cout << "> accumulators: " << std::endl << accumulators << std::endl;

	{
		alps::hdf5::archive ar("test.h5", "w");
		ar["/measurements"] << accumulators;
	}

	{
		alps::hdf5::archive ar("test.h5", "r");
		accumulator_set accumulators2;
		ar["/measurements"] >> accumulators2;

	    BOOST_REQUIRE(count(accumulators2["mean"]) == 5);
	    BOOST_REQUIRE(mean(accumulators2["mean"].get<double>()) == 6.2);

		std::cout << "> accumulators: " << std::endl;
	    for (accumulator_set::const_iterator it = accumulators.begin(); it != accumulators.end(); ++it) {
	    	std::stringstream sa;
	    	sa << *it->second;
	    	std::stringstream sa2;
	    	sa2 << accumulators2[it->first];

	    	std::cout << it->first << ": " << sa.str() << " <=> " << sa2.str() << std::endl;
		    BOOST_REQUIRE(sa.str() == sa2.str());
	    }
	    std::cout << std::endl;

	    accumulators2.reset();
	    accumulators2["mean"] << 1;

	    BOOST_REQUIRE(count(accumulators2["mean"]) == 1);
	    BOOST_REQUIRE(mean(accumulators2["mean"].get<double>()) == 1);
	}

	result_set results(accumulators);

    BOOST_REQUIRE(count(results["mean"]) == 5);
    BOOST_REQUIRE(mean(results["mean"].get<double>()) == 6.2);

	std::cout << "> results: " << std::endl << results << std::endl;

	std::cout << "> sin(\"vectors\"): " << std::endl << sin(results["vector"]) << std::endl << std::endl;

	{
		alps::hdf5::archive ar("test.h5", "w");
		ar["/results"] << results;
	}

	{
		alps::hdf5::archive ar("test.h5", "r");
		accumulator_set results2;
		ar["/results"] >> results2;

	    BOOST_REQUIRE(count(results2["mean"]) == 5);
	    BOOST_REQUIRE(mean(results2["mean"].get<double>()) == 6.2);

		std::cout << "> results: " << std::endl;
	    for (result_set::const_iterator it = results.begin(); it != results.end(); ++it) {
	    	std::stringstream sa;
	    	sa << *it->second;
	    	std::stringstream sa2;
	    	sa2 << results2[it->first];

	    	std::cout << it->first << ": " << sa.str() << " <=> " << sa2.str() << std::endl;
		    BOOST_REQUIRE(sa.str() == sa2.str());
	    }
	    std::cout << std::endl;

	}

/* TODO:
- implement functions for results
- implement jacknife for results
- implement boost::ArgPack for external weight
- implement tag-magic
- add has hdf5
- plug it into the easy examples
*/
}
