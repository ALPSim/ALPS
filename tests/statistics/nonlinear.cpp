// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Error of a nonlinear function of two correlated averages. For x ~ N(0, 1),
// A = x^2 and B = x^4, the ratio U = <B> / <A>^2 is 3, and the delta method
// with Var(A) = 2, Var(B) = 96, Cov(A, B) = 12 gives Var(U) = 24 / N.
// Ignoring the covariance would give 168 / N, an error 2.65 times too large.
#include "streams.hpp"

#include <alps/alea.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <vector>

namespace {

const std::size_t n = std::size_t(1) << 16;
const double exact = 3.0;
const double exact_error = std::sqrt(24.0 / n);

// A jackknife over 128 bins scatters by 6%, and the high moments of x^4 add to
// that: 25% is about four standard deviations.
void expect_ratio(double value, double error, const char * library) {
    EXPECT_NEAR(exact, value, 4 * exact_error) << library;
    EXPECT_NEAR(1.0, error / exact_error, 0.25) << library;
}

double ratio(double a, double b) { return b / (a * a); }

}

TEST(Nonlinear, Alea)
{
    using namespace alps::alea;
    batch_acc<double> batch(2);
    cov_acc<double> cov(2);
    var_acc<double> var(2);
    statistics_test::normal xi(424242);
    for (std::size_t t = 0; t < n; ++t) {
        const double x = xi();
        const std::vector<double> ab = {x * x, x * x * x * x};
        batch << ab;
        cov << ab;
        var << ab;
    }
    const auto tf = make_transformer(std::function<double(double, double)>(ratio));

    const batch_result<double> jackknife = transform(jackknife_prop(), tf, batch.finalize());
    expect_ratio(jackknife.mean()[0], jackknife.stderror()[0], "ALEA v2 jackknife");

    const cov_result<double> linear = transform(linear_prop(), tf, cov.finalize());
    expect_ratio(linear.mean()[0], linear.stderror()[0], "ALEA v2 linear propagation");

    // Without the covariance the same propagation is 2.65 times too large;
    // this is why results must not be combined as if they were independent.
    const cov_result<double> uncorrelated = transform(linear_prop(), tf, var.finalize());
    EXPECT_NEAR(std::sqrt(168.0 / 24.0), uncorrelated.stderror()[0] / exact_error, 0.4);
}
