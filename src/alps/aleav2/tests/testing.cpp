/*
 * Copyright (C) 2026 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */
// Hotelling's T^2 test (testing.hpp) against hand-computed statistics, and
// its false-alarm rate on data where the null hypothesis holds.
#include <alps/alea.hpp>
#include <alps/alea/testing.hpp>

#include "gtest/gtest.h"

#include <Eigen/Dense>
#include <cmath>
#include <random>

namespace {

using alps::alea::cov_acc;
using alps::alea::cov_result;

// Correlated Gaussian vectors with covariance l * l^T around `mean`
cov_result<double> sample(std::mt19937 &engine, const Eigen::VectorXd &mean,
                          const Eigen::MatrixXd &l, size_t count)
{
    std::normal_distribution<double> normal;
    cov_acc<double> acc(mean.size());
    for (size_t i = 0; i != count; ++i) {
        Eigen::VectorXd xi(mean.size());
        for (int k = 0; k != xi.size(); ++k)
            xi[k] = normal(engine);
        acc << Eigen::VectorXd(mean + l * xi);
    }
    return acc.finalize();
}

Eigen::MatrixXd mixing(int n, double scale)
{
    // Lower triangular, so that the components are correlated and the
    // principal axes are not the coordinate axes.
    Eigen::MatrixXd l = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i != n; ++i)
        for (int j = 0; j <= i; ++j)
            l(i, j) = scale * (i == j ? 1.0 : 0.5 / (1 + i - j));
    return l;
}

} // namespace

// One sample with a full covariance: T^2 = N d^T C^-1 d, referred to
// F(n, N - n) after scaling by (N - n) / (n (N - 1)) (Wallerberger and Gull,
// Eq. 9). The rotation into the principal axes must not change it.
TEST(t2_test, one_sample_covariance_matches_hotelling)
{
    std::mt19937 engine(7);
    Eigen::VectorXd mean(3);
    mean << 1.0, -2.0, 0.5;
    cov_result<double> res = sample(engine, mean, mixing(3, 1.0), 20);

    Eigen::VectorXd expected = mean;
    expected[1] += 0.3;
    const double n = 3, nobs = res.observations();
    const Eigen::VectorXd d = res.mean() - expected;
    const double t2 = nobs * d.dot(res.cov().ldlt().solve(d));

    alps::alea::t2_result r = alps::alea::test_mean(res, expected);
    EXPECT_NEAR(r.score(), (nobs - n) / (n * (nobs - 1)) * t2, 1e-10 * r.score());
    EXPECT_EQ(r.dist().degrees_of_freedom1(), n);
    EXPECT_EQ(r.dist().degrees_of_freedom2(), nobs - n);
}

// Two samples of different length and covariance: Welch form S = C1/N1 +
// C2/N2 with the Krishnamoorthy-Yu degrees of freedom, written here with
// traces of S^-1 A_k instead of the eigenbasis the library uses.
TEST(t2_test, two_sample_welch_matches_krishnamoorthy_yu)
{
    std::mt19937 engine(11);
    Eigen::VectorXd mean(3);
    mean << 1.0, -2.0, 0.5;
    cov_result<double> r1 = sample(engine, mean, mixing(3, 3.0), 8);
    cov_result<double> r2 = sample(engine, mean, mixing(3, 0.5), 32);

    const double n = 3, n1 = r1.observations(), n2 = r2.observations();
    const Eigen::MatrixXd a1 = r1.cov() / n1, a2 = r2.cov() / n2;
    const Eigen::MatrixXd s = a1 + a2;
    const Eigen::VectorXd d = r1.mean() - r2.mean();
    const double t2 = d.dot(s.ldlt().solve(d));

    const Eigen::MatrixXd b1 = s.ldlt().solve(a1), b2 = s.ldlt().solve(a2);
    const double denom = ((b1 * b1).trace() + b1.trace() * b1.trace()) / (n1 - 1)
                       + ((b2 * b2).trace() + b2.trace() * b2.trace()) / (n2 - 1);
    const double nu = (n + n * n) / denom;

    alps::alea::t2_result r = alps::alea::test_mean(r1, r2);
    EXPECT_NEAR(r.dist().degrees_of_freedom2(), nu + 1 - n, 1e-10 * nu);
    EXPECT_NEAR(r.score(), (nu + 1 - n) / (n * nu) * t2, 1e-10 * r.score());

    // A shift of many standard errors is rejected
    Eigen::VectorXd shifted = mean;
    shifted[0] += 30.0;  // about 28 standard errors of r1
    cov_result<double> r3 = sample(engine, shifted, mixing(3, 0.5), 32);
    EXPECT_LT(alps::alea::test_mean(r1, r3).pvalue(), 1e-4);
}

// One component: the Welch t test, with Welch-Satterthwaite degrees of freedom
TEST(t2_test, two_sample_scalar_is_welch_t)
{
    std::mt19937 engine(3);
    Eigen::VectorXd mean(1);
    mean << 2.0;
    Eigen::MatrixXd wide(1, 1), narrow(1, 1);
    wide << 4.0;
    narrow << 1.0;
    cov_result<double> r1 = sample(engine, mean, wide, 10);
    cov_result<double> r2 = sample(engine, mean, narrow, 40);

    const double a1 = r1.cov()(0, 0) / r1.observations();
    const double a2 = r2.cov()(0, 0) / r2.observations();
    const double nu = (a1 + a2) * (a1 + a2)
                    / (a1 * a1 / (r1.observations() - 1) + a2 * a2 / (r2.observations() - 1));
    const double t = (r1.mean()[0] - r2.mean()[0]) / std::sqrt(a1 + a2);

    alps::alea::t2_result r = alps::alea::test_mean(r1, r2);
    EXPECT_NEAR(r.dist().degrees_of_freedom2(), nu, 1e-10 * nu);
    EXPECT_NEAR(r.score(), t * t, 1e-10 * t * t);
}

TEST(t2_result, pvalue_is_two_sided)
{
    // Lower tail testable (more than three components)
    alps::alea::t2_result low(0.2, 5, 20), high(3.0, 5, 20), mid(0.9, 5, 20);
    EXPECT_DOUBLE_EQ(low.pvalue(), 2 * low.pvalue_lower());
    EXPECT_DOUBLE_EQ(high.pvalue(), 2 * high.pvalue_upper());
    EXPECT_DOUBLE_EQ(mid.pvalue(), std::min(1.0, 2 * std::min(mid.pvalue_lower(),
                                                              mid.pvalue_upper())));
    // Lower tail not testable: one-sided upper test
    alps::alea::t2_result few(3.0, 2, 10);
    EXPECT_DOUBLE_EQ(few.pvalue(), few.pvalue_upper());
    // Too few observations: NaN, so that a pvalue() >= alpha check fails
    alps::alea::t2_result none(NAN, 2, 0.5);
    EXPECT_TRUE(std::isnan(none.pvalue()));
}

// Exactly correlated components make the covariance singular: here
// (x, y, x, y, x + y, x - y), rank 2, with a mean of 30 standard deviations,
// so that the null eigenvalues are rounding of either sign. They must be
// dropped whatever the units of the data, for one and for two samples.
cov_result<double> sum_rule_sample(std::mt19937 &engine, double scale, int count)
{
    std::normal_distribution<double> normal;
    cov_acc<double> acc(6);
    for (int i = 0; i != count; ++i) {
        const double x = scale * (30 + normal(engine)), y = scale * (-20 + normal(engine));
        Eigen::VectorXd v(6);
        v << x, y, x, y, x + y, x - y;
        acc << v;
    }
    return acc.finalize();
}

TEST(t2_test, exactly_correlated_components_are_dropped)
{
    Eigen::VectorXd exact(6);
    exact << 30, -20, 30, -20, 10, 50;
    for (double scale : {1e-9, 1.0, 1e9}) {
        std::mt19937 engine(5);
        for (int rep = 0; rep != 200; ++rep) {
            cov_result<double> r1 = sum_rule_sample(engine, scale, 50);
            cov_result<double> r2 = sum_rule_sample(engine, scale, 80);
            alps::alea::t2_result one = alps::alea::test_mean(r1, Eigen::VectorXd(scale * exact));
            alps::alea::t2_result two = alps::alea::test_mean(r1, r2);
            ASSERT_EQ(one.dist().degrees_of_freedom1(), 2) << "scale " << scale << " rep " << rep;
            ASSERT_EQ(two.dist().degrees_of_freedom1(), 2) << "scale " << scale << " rep " << rep;
            ASSERT_TRUE(std::isfinite(one.score()) && one.score() >= 0) << "scale " << scale;
            ASSERT_TRUE(std::isfinite(two.score()) && two.score() >= 0) << "scale " << scale;
        }
    }
}

// False-alarm rate where the null hypothesis holds exactly (independent
// Gaussian samples): 2000 tests at alpha = 0.05 reject 100 +- 9.7. The band
// [51, 149] is five standard deviations; a correct test fails it with
// probability below 1e-6. Taking the smaller tail without doubling rejects
// about 200; pooling unequal variances rejects far more than 100.
TEST(t2_test, false_alarm_rate_one_sample)
{
    std::mt19937 engine(42);
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(5);
    const Eigen::MatrixXd l = mixing(5, 1.0);
    int rejected = 0;
    for (int rep = 0; rep != 2000; ++rep) {
        cov_result<double> res = sample(engine, mean, l, 30);
        if (alps::alea::test_mean(res, mean).pvalue() < 0.05)
            ++rejected;
    }
    EXPECT_GE(rejected, 51);
    EXPECT_LE(rejected, 149);
}

TEST(t2_test, false_alarm_rate_two_sample_unequal_variance)
{
    std::mt19937 engine(43);
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(4);
    const Eigen::MatrixXd wide = mixing(4, 4.0), narrow = mixing(4, 1.0);
    int rejected = 0;
    for (int rep = 0; rep != 2000; ++rep) {
        cov_result<double> r1 = sample(engine, mean, wide, 10);
        cov_result<double> r2 = sample(engine, mean, narrow, 40);
        if (alps::alea::test_mean(r1, r2).pvalue() < 0.05)
            ++rejected;
    }
    EXPECT_GE(rejected, 51);
    EXPECT_LE(rejected, 149);
}
