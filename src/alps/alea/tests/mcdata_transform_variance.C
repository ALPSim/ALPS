/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* Copyright (C) 1994-2025 by the ALPS collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Regression test for the variance_opt self-assignment bug in mcdata::transform_linear
// and mcdata::transform. Each method accepted a variance_opt parameter but assigned
// variance_opt_ to itself instead of to the parameter.

#include <alps/alea.h>
#include <alps/alea/mcdata.hpp>
#include <boost/optional.hpp>
#include <gtest/gtest.h>

class McdataVariance : public ::testing::Test {
protected:
    alps::RealObservable observable{"test"};
    void SetUp() override {
        for (int i = 0; i < 1024; ++i) observable << double(i % 10);
    }
};

TEST_F(McdataVariance, LinearTransformStoresSuppliedVariance)
{
    alps::alea::mcdata<double> data(observable);
    data.transform_linear([](double x) { return 2. * x; }, 1.5,
                          boost::optional<double>(42.));
    ASSERT_TRUE(data.has_variance());
    EXPECT_NEAR(data.variance(), 42., 1e-10);
}

TEST_F(McdataVariance, LinearTransformCanClearVariance)
{
    alps::alea::mcdata<double> data(observable);
    data.transform_linear([](double x) { return x; }, 1.5, boost::optional<double>(42.));
    data.transform_linear([](double x) { return x; }, 1.5, boost::none);
    EXPECT_FALSE(data.has_variance());
}

TEST_F(McdataVariance, UnaryTransformStoresSuppliedVariance)
{
    alps::alea::mcdata<double> data(observable);
    data.transform([](double x) { return x + 1.; }, 2., boost::optional<double>(99.));
    ASSERT_TRUE(data.has_variance());
    EXPECT_NEAR(data.variance(), 99., 1e-10);
}

TEST_F(McdataVariance, BinaryTransformStoresSuppliedVariance)
{
    alps::alea::mcdata<double> lhs(observable), rhs(observable);
    lhs.transform(rhs, [](double x, double y) { return x + y; }, 3.,
                  boost::optional<double>(7.));
    ASSERT_TRUE(lhs.has_variance());
    EXPECT_NEAR(lhs.variance(), 7., 1e-10);
}
