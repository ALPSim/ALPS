// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Coverage of the reported error bars: over many independent streams with
// true mean 0, the fraction of streams whose mean lies within k reported
// errors must be the Gaussian value, 0.6827 for k = 1 and 0.9545 for k = 2.
// The mean of a Gaussian AR(1) series is exactly Gaussian. The band is 3.5
// binomial standard deviations; seeds are fixed, so the outcome is too.
#include "streams.hpp"

#include <alps/alea.hpp>

#include <gtest/gtest.h>

#include <cmath>

namespace {

const int streams = 500;
const std::size_t n = std::size_t(1) << 15;

struct counts {
    int within[2] = {0, 0};
    void add(double mean, double error) {
        within[0] += std::abs(mean) < error;
        within[1] += std::abs(mean) < 2 * error;
    }
    double fraction(int k) const { return double(within[k - 1]) / streams; }
};

struct coverage {
    counts batch, autocorr;
};

coverage run(double phi) {
    coverage result;
    for (int r = 0; r < streams; ++r) {
        alps::alea::batch_acc<double> batch(1);
        alps::alea::autocorr_acc<double> autocorr(1);
        statistics_test::ar1 source(7919u * r + 1, phi);
        for (std::size_t t = 0; t < n; ++t) {
            const double x = source();
            batch << x;
            autocorr << x;
        }
        const alps::alea::batch_result<double> batch_result = batch.finalize();
        result.batch.add(batch_result.mean()[0], batch_result.stderror()[0]);
        const alps::alea::autocorr_result<double> autocorr_result = autocorr.finalize();
        result.autocorr.add(autocorr_result.mean()[0], autocorr_result.stderror()[0]);
    }
    return result;
}

const double nominal[2] = {0.6827, 0.9545};

double band(int k) { return 3.5 * std::sqrt(nominal[k - 1] * (1 - nominal[k - 1]) / streams); }

void expect_nominal(const counts & c, const char * library) {
    for (int k = 1; k <= 2; ++k)
        EXPECT_NEAR(nominal[k - 1], c.fraction(k), band(k)) << library << ", k = " << k;
}

}

TEST(Coverage, Uncorrelated)
{
    const coverage result = run(0.0);
    expect_nominal(result.batch, "ALEA v2 batch");
    expect_nominal(result.autocorr, "ALEA v2 autocorr");
}

TEST(Coverage, Correlated)
{
    const coverage result = run(0.9);
    expect_nominal(result.batch, "ALEA v2 batch");
    expect_nominal(result.autocorr, "ALEA v2 autocorr");
}
