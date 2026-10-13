/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: detailedbinning.C 3654 2010-01-06 23:47:46Z troyer $ */

#include <alps/alea.h>
#include <alps/alea/abstractsimpleobservable.ipp>
#include <alps/alea/simpleobservable.ipp>
#include <alps/alea/simpleobseval.ipp>
#include <boost/random.hpp>
#include <gtest/gtest.h>
#include <complex>

TEST(ComplexObservable, MeansMatchIndependentSums)
{
    boost::minstd_rand0 engine(1u);
    boost::uniform_01<boost::minstd_rand0> random(engine);
    alps::ComplexObservable a("Complex Observable A"), b("Complex Observable B");
    std::complex<double> sum_a{}, sum_b{};
    for (int i = 0; i < 1000; ++i) {
        const double ar = random(), ai = random() + 1.;
        const double br = random() + 3., bi = random() + 1.;
        const std::complex<double> av(ar, ai), bv(br, bi);
        a << av;
        b << bv;
        sum_a += av;
        sum_b += bv;
    }
    EXPECT_EQ(a.count(), 1000u);
    EXPECT_EQ(b.count(), 1000u);
    EXPECT_NEAR(std::abs(a.mean() - sum_a / 1000.), 0., 1e-12);
    EXPECT_NEAR(std::abs(b.mean() - sum_b / 1000.), 0., 1e-12);
    const alps::ComplexObsevaluator eval_a(a), eval_b(b);
    const alps::ComplexObsevaluator product = eval_a * eval_b;
    EXPECT_EQ(product.count(), 1000u);
    EXPECT_TRUE(std::isfinite(product.mean().real()));
    EXPECT_TRUE(std::isfinite(product.mean().imag()));
    EXPECT_TRUE(std::isfinite(product.error().real()));
    EXPECT_TRUE(std::isfinite(product.error().imag()));
}

TEST(ComplexObservable, ConstantProductHasAnalyticMeanAndZeroError)
{
    alps::ComplexObservable a("A"), b("B");
    for (int i = 0; i < 1024; ++i) {
        a << std::complex<double>(1., 2.);
        b << std::complex<double>(3., 1.);
    }
    const alps::ComplexObsevaluator eval_a(a), eval_b(b);
    const alps::ComplexObsevaluator product = eval_a * eval_b;
    EXPECT_NEAR(std::abs(product.mean() - std::complex<double>(1., 7.)), 0., 1e-12);
    EXPECT_NEAR(std::abs(product.error()), 0., 1e-12);
}
