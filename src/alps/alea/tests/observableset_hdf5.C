/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010-2012 by Lukas Gamper <gamperl -at- gmail.com>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/alea.h>
#include <alps/hdf5.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <boost/random.hpp>
#include "observable_checks.hpp"
#include <cmath>

namespace {
void expect_measurements(const alps::ObservableSet& measurement)
{
    EXPECT_EQ(measurement.get<alps::RealObservable>("No Measurements").count(), 0u);
    const auto& sign = measurement.get<alps::RealObservable>("Sign");
    const auto& second = measurement.get<alps::RealObservable>("Test 2");
    const auto& third = measurement.get<alps::RealObservable>("Test 3");
    const auto& fourth = measurement.get<alps::SimpleRealObservable>("Test 4");
    const auto& signed_observable = measurement.get<alps::SignedObservable<alps::RealObservable>>("Test");
    const auto& raw = dynamic_cast<const alps::RealObservable&>(signed_observable.signed_observable());
    const alps::RealObsevaluator ratio = measurement["Ratio"];
    alps_test::expect_estimate(sign, {{1., 0., 0., 0.}}, "Sign");
    alps_test::expect_estimate(second, {{0.493208, 0.00277, 5.00001e-7, 5.00001e-6}}, "Test 2");
    alps_test::expect_estimate(third, {{0.507869, 0.00321, 5.00001e-7, 5.00001e-6}}, "Test 3");
    alps_test::expect_estimate(fourth, {{0.498397, 0.00287, 5.00001e-7, 5.00001e-6}}, "Test 4");
    alps_test::expect_estimate(signed_observable, {{0.500296, 0.0030556, 5.00001e-7, 5.00001e-8}}, "Test");
    alps_test::expect_estimate(raw, {{0.500301, 0.00282, 5.00001e-7, 5.00001e-6}}, "Sign * Test");
    alps_test::expect_estimate(ratio, {{0.970608, 0.00847, 5.00001e-7, 5.00001e-6}}, "Ratio");
    EXPECT_EQ(signed_observable.sign_name(), "Sign");
    EXPECT_EQ(third.converged_errors(), alps::MAYBE_CONVERGED);
    // The former formatter printed zero for this undefined constant-data tau.
    EXPECT_TRUE(std::isnan(sign.tau()));
    EXPECT_NEAR(raw.tau(), -0.0227, 5.00001e-5);
    EXPECT_NEAR(second.tau(), -0.0374, 5.00001e-5);
    EXPECT_NEAR(third.tau(), 0.116, 5.00001e-4);
    const double raw_errors[] = {0.00288344, 0.00287294, 0.00288716, 0.00288336, 0.00291027, 0.00294589, 0.00281725};
    const double second_errors[] = {0.0028763, 0.00288887, 0.00290494, 0.00287583, 0.00280704, 0.00282259, 0.00276664};
    const double third_errors[] = {0.00289105, 0.00286199, 0.00287098, 0.00289687, 0.00288711, 0.00288283, 0.00320908};
    EXPECT_EQ(sign.count(), 10000u);
    EXPECT_EQ(raw.count(), 10000u);
    EXPECT_EQ(second.count(), 10000u);
    EXPECT_EQ(third.count(), 10000u);
    EXPECT_EQ(fourth.count(), 10000u);
    for (unsigned level = 0; level < 7; ++level) {
        SCOPED_TRACE(::testing::Message() << "bin level " << level);
        EXPECT_EQ(sign.error(level), 0.);
        EXPECT_NEAR(raw.error(level), raw_errors[level], 5.00001e-9);
        EXPECT_NEAR(second.error(level), second_errors[level], 5.00001e-9);
        EXPECT_NEAR(third.error(level), third_errors[level], 5.00001e-9);
    }
    const alps::IntHistogramObsevaluator histogram = measurement["Histogram"];
    const unsigned expected[] = {966, 1012, 998, 988, 951, 1060, 1060, 988, 1009, 968};
    EXPECT_EQ(histogram.count(), 10000u);
    ASSERT_EQ(histogram.size(), 10u);
    for (int bin = 0; bin < 10; ++bin) EXPECT_EQ(histogram[bin], expected[bin]);
}
}

TEST(ObservableSetHdf5, RestoresEmptyAndPopulatedSets)
{
    alps::testing::TemporaryDirectory temporary;
    const auto filename = (temporary.path() / "observableset.h5").string();
    {
        alps::ObservableSet measurement;
        measurement << alps::make_observable(alps::RealObservable("Test"), true)
                    << alps::RealObservable("Sign")
                    << alps::RealObservable("No Measurements")
                    << alps::IntHistogramObservable("Histogram", 0, 10)
                    << alps::RealObservable("Test 2")
                    << alps::RealObservable("Test 3")
                    << alps::SimpleRealObservable("Test 4");
        alps::hdf5::archive archive(filename, "a");
        archive["/test/0/result"] << measurement;
    }
    {
        boost::minstd_rand0 engine(1u);
        boost::uniform_01<boost::minstd_rand0> random(engine);
        alps::ObservableSet measurement;
        {
            alps::hdf5::archive archive(filename, "r");
            archive["/test/0/result"] >> measurement;
        }
        EXPECT_EQ(measurement.get<alps::RealObservable>("Sign").count(), 0u);
        for (int i = 0; i < 10000; ++i) {
            measurement["Test"] << random();
            measurement["Sign"] << 1.;
            measurement["Histogram"] << static_cast<int>(10 * random());
            measurement["Test 2"] << random();
            measurement["Test 3"] << random();
            measurement["Test 4"] << random();
        }
        const alps::RealObsevaluator second = measurement["Test 2"], third = measurement["Test 3"];
        alps::RealObsevaluator ratio("Ratio");
        ratio = second / third;
        measurement.addObservable(ratio);
        alps::hdf5::archive archive(filename, "a");
        archive["/test/0/result"] << measurement;
    }
    for (bool prepopulate : {false, true}) {
        SCOPED_TRACE(::testing::Message() << "prepopulate " << prepopulate);
        alps::ObservableSet measurement;
        if (prepopulate)
            measurement << alps::make_observable(alps::RealObservable("Test"), true)
                        << alps::RealObservable("Sign")
                        << alps::RealObservable("No Measurements")
                        << alps::IntHistogramObservable("Histogram", 0, 10);
        alps::hdf5::archive archive(filename, "r");
        archive["/test/0/result"] >> measurement;
        expect_measurements(measurement);
    }
}
