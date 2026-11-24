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
#include <functional>
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

// The cut for null directions is taken per component, on the covariance
// scaled by D = diag(sqrt(var + mean^2)). Independent components with very
// different scales are all tested: here a component of mean 1000 and
// standard deviation 1 next to one of mean 0 and standard deviation `small`.
// 1000 tests at alpha = 0.05 reject 50 +- 6.9; [16, 84] is five standard
// deviations. A cut taken from the largest component treated the small one
// as null and rejected every test.
TEST(t2_test, components_of_different_scale_are_tested)
{
    for (double small : {1e-2, 1e-5, 1e-8, 1e-12}) {
        std::mt19937 engine(3);
        Eigen::VectorXd mean(2);
        mean << 1000.0, 0.0;
        Eigen::MatrixXd l = Eigen::MatrixXd::Zero(2, 2);
        l(0, 0) = 1.0;
        l(1, 1) = small;
        int rejected_one = 0, rejected_two = 0;
        for (int rep = 0; rep != 1000; ++rep) {
            cov_result<double> r1 = sample(engine, mean, l, 100);
            cov_result<double> r2 = sample(engine, mean, l, 300);
            alps::alea::t2_result one = alps::alea::test_mean(r1, mean);
            alps::alea::t2_result two = alps::alea::test_mean(r1, r2);
            ASSERT_EQ(one.dist().degrees_of_freedom1(), 2) << "std " << small;
            ASSERT_EQ(two.dist().degrees_of_freedom1(), 2) << "std " << small;
            rejected_one += one.pvalue() < 0.05;
            rejected_two += two.pvalue() < 0.05;
        }
        EXPECT_GE(rejected_one, 16) << "std " << small;
        EXPECT_LE(rejected_one, 84) << "std " << small;
        EXPECT_GE(rejected_two, 16) << "std " << small;
        EXPECT_LE(rejected_two, 84) << "std " << small;
    }
}

// A constant component that is not exactly representable (0.1) is not
// exact after accumulation: the error of its mean and variance grows like
// count * eps (not sqrt(count) * eps, since every addition rounds the same
// way). It must be dropped as exact, for a variance and for a covariance,
// and the remaining component tested as usual: 25 tests at alpha = 0.05
// reject 1.25 +- 1.1. The data cannot tell such a constant from one that
// fluctuates below the precision of its variance (here 7e-6 relative), so
// a constant is resolved to about 2e-8 relative: one off by 1e-6 must be
// rejected.
TEST(t2_test, constant_components_are_dropped_after_rounding)
{
    using alps::alea::var_acc;
    std::mt19937 engine(5);
    std::normal_distribution<double> normal;
    Eigen::VectorXd expected(2), wrong(2);
    expected << 0.1, 0.0;
    wrong << 0.1 * (1 + 1e-6), 0.0;
    int rejected_var = 0, rejected_cov = 0, rejected_two = 0;
    for (int rep = 0; rep != 25; ++rep) {
        var_acc<double> vacc(2), vacc2(2);
        cov_acc<double> cacc(2);
        for (int i = 0; i != 100000; ++i) {
            Eigen::VectorXd x(2);
            x << 0.1, normal(engine);
            vacc << x;
            cacc << x;
            x[1] = normal(engine);
            vacc2 << x;
        }
        alps::alea::var_result<double> vres = vacc.finalize(), vres2 = vacc2.finalize();
        cov_result<double> cres = cacc.finalize();

        alps::alea::t2_result tv = alps::alea::test_mean(vres, expected);
        alps::alea::t2_result tc = alps::alea::test_mean(cres, expected);
        alps::alea::t2_result t2 = alps::alea::test_mean(vres, vres2);
        ASSERT_EQ(tv.dist().degrees_of_freedom1(), 1) << "rep " << rep;
        ASSERT_EQ(tc.dist().degrees_of_freedom1(), 1) << "rep " << rep;
        ASSERT_EQ(t2.dist().degrees_of_freedom1(), 1) << "rep " << rep;
        rejected_var += tv.pvalue() < 0.05;
        rejected_cov += tc.pvalue() < 0.05;
        rejected_two += t2.pvalue() < 0.05;

        if (rep == 0) {
            EXPECT_LT(alps::alea::test_mean(vres, wrong).pvalue(), 1e-10);
            EXPECT_LT(alps::alea::test_mean(cres, wrong).pvalue(), 1e-10);
        }
    }
    EXPECT_LE(rejected_var, 8);
    EXPECT_LE(rejected_cov, 8);
    EXPECT_LE(rejected_two, 8);
}

// Results without any variance: nothing is left to test if they agree, so
// the test passes (p = 1); if they differ, it rejects.
TEST(t2_test, exact_results)
{
    Eigen::VectorXd x(2), y(2);
    x << 1.0, 2.0;
    y << 1.0, 2.5;
    cov_acc<double> a(2), b(2), c(2);
    alps::alea::var_acc<double> va(2);
    for (int i = 0; i != 10; ++i) {
        a << x;
        b << x;
        c << y;
        va << x;
    }
    cov_result<double> ra = a.finalize(), rb = b.finalize(), rc = c.finalize();
    alps::alea::var_result<double> rva = va.finalize();

    EXPECT_EQ(alps::alea::test_mean(ra, x).pvalue(), 1.0);
    EXPECT_EQ(alps::alea::test_mean(rva, x).pvalue(), 1.0);
    EXPECT_EQ(alps::alea::test_mean(ra, rb).pvalue(), 1.0);
    EXPECT_EQ(alps::alea::test_mean(ra, rva).pvalue(), 1.0);
    EXPECT_LT(alps::alea::test_mean(ra, y).pvalue(), 1e-10);
    EXPECT_LT(alps::alea::test_mean(rva, y).pvalue(), 1e-10);
    EXPECT_LT(alps::alea::test_mean(ra, rc).pvalue(), 1e-10);
}

// A NaN in the data must give a NaN p-value (so that `pvalue() >= alpha`
// fails), not be dropped. Empty results have nothing to test.
TEST(t2_test, nan_propagates_and_empty_passes)
{
    std::mt19937 engine(9);
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(3);
    cov_result<double> good = sample(engine, mean, mixing(3, 1.0), 50);

    cov_acc<double> acc(3);
    alps::alea::var_acc<double> vacc(3);
    std::normal_distribution<double> normal;
    for (int i = 0; i != 50; ++i) {
        Eigen::VectorXd v(3);
        v << normal(engine), normal(engine), i == 7 ? NAN : normal(engine);
        acc << v;
        vacc << v;
    }
    cov_result<double> bad = acc.finalize();
    alps::alea::var_result<double> vbad = vacc.finalize();
    EXPECT_TRUE(std::isnan(alps::alea::test_mean(bad, mean).pvalue()));
    EXPECT_TRUE(std::isnan(alps::alea::test_mean(vbad, mean).pvalue()));
    EXPECT_TRUE(std::isnan(alps::alea::test_mean(bad, good).pvalue()));
    EXPECT_TRUE(std::isnan(alps::alea::test_mean(vbad, good).pvalue()));

    Eigen::VectorXd nan_expected = mean;
    nan_expected[1] = NAN;
    EXPECT_TRUE(std::isnan(alps::alea::test_mean(good, nan_expected).pvalue()));

    cov_acc<double> e1(0), e2(0);
    alps::alea::var_acc<double> e3(0);
    for (int i = 0; i != 5; ++i) {
        e1 << Eigen::VectorXd(0);
        e2 << Eigen::VectorXd(0);
        e3 << Eigen::VectorXd(0);
    }
    cov_result<double> re1 = e1.finalize(), re2 = e2.finalize();
    alps::alea::var_result<double> re3 = e3.finalize();
    EXPECT_EQ(alps::alea::test_mean(re1, Eigen::VectorXd(0)).pvalue(), 1.0);
    EXPECT_EQ(alps::alea::test_mean(re3, Eigen::VectorXd(0)).pvalue(), 1.0);
    EXPECT_EQ(alps::alea::test_mean(re1, re2).pvalue(), 1.0);
}

// std::min(1.0, x) is 1 for x = NaN: pvalue() must return NaN itself when the
// two-sided p-value is taken (more than three components).
TEST(t2_result, nan_score_gives_nan_pvalue_for_many_components)
{
    EXPECT_TRUE(std::isnan(alps::alea::t2_result(NAN, 5, 20).pvalue()));
    EXPECT_TRUE(std::isnan(alps::alea::t2_result(NAN, 2, 20).pvalue()));
}

// rtol is the relative accuracy of exact components. Here the second
// component equals the first by symmetry, but the measurement computes it
// only to about 1e-6 relative (a bias of 1e-6 and noise of 1e-7). With the
// default rtol (64 eps) the broken symmetry is resolved and the test rejects;
// with rtol = 1e-5 the components are exact to the stated accuracy, the
// direction is dropped, and the remaining one is tested as usual.
TEST(t2_test, rtol_is_the_accuracy_of_exact_components)
{
    std::mt19937 engine(17);
    std::normal_distribution<double> normal;
    cov_acc<double> acc(2);
    for (int i = 0; i != 10000; ++i) {
        const double x = 3 + normal(engine);
        Eigen::VectorXd v(2);
        v << x, x * (1 + 1e-6) + 1e-7 * normal(engine);
        acc << v;
    }
    cov_result<double> res = acc.finalize();
    Eigen::VectorXd expected(2);
    expected << 3.0, 3.0;
    EXPECT_LT(alps::alea::test_mean(res, expected).pvalue(), 1e-6);
    alps::alea::t2_result loose = alps::alea::test_mean(res, expected, 1e-5);
    EXPECT_EQ(loose.dist().degrees_of_freedom1(), 1);
    EXPECT_GE(loose.pvalue(), 0.01);
}

// atol is in the units of the data (not squared): a component that should be
// zero and is computed as 1e-12 is exact for atol = 1e-11 and not for 1e-13.
TEST(t2_test, atol_is_in_units_of_the_data)
{
    std::mt19937 engine(19);
    std::normal_distribution<double> normal;
    alps::alea::var_acc<double> acc(2);
    for (int i = 0; i != 1000; ++i) {
        Eigen::VectorXd v(2);
        v << normal(engine), 1e-12;
        acc << v;
    }
    alps::alea::var_result<double> res = acc.finalize();
    const Eigen::VectorXd zero = Eigen::VectorXd::Zero(2);
    alps::alea::t2_result loose = alps::alea::test_mean(res, zero, 0, 1e-11);
    EXPECT_EQ(loose.dist().degrees_of_freedom1(), 1);
    EXPECT_GE(loose.pvalue(), 1e-4);
    EXPECT_LT(alps::alea::test_mean(res, zero, 0, 1e-13).pvalue(), 1e-6);
}

// Fewer observations than components: Hotelling's T^2 does not exist, and
// the test of Srivastava and Du is used. Its false-alarm rate at alpha = 0.05
// is approximate (0.03 to 0.08 in simulations with at least 15 observations),
// so the band here is [20, 100] of 1000, 2% to 10%. Before, the result was NaN
// for one sample. For two samples, the smaller one here has fewer
// observations than components and the larger variance: Hotelling's T^2 with
// the Krishnamoorthy-Yu degrees of freedom rejected 47% of these tests.
namespace {

int rejected_of_1000(const std::function<alps::alea::t2_result(std::mt19937 &)> &test)
{
    std::mt19937 engine(23);
    int rejected = 0;
    for (int rep = 0; rep != 1000; ++rep) {
        alps::alea::t2_result r = test(engine);
        EXPECT_TRUE(r.normal());
        EXPECT_FALSE(std::isnan(r.pvalue()));
        rejected += r.pvalue() < 0.05;
    }
    return rejected;
}

alps::alea::var_result<double> independent_sample(std::mt19937 &engine, int p, double sd, size_t count)
{
    std::normal_distribution<double> normal;
    alps::alea::var_acc<double> acc(p);
    for (size_t i = 0; i != count; ++i) {
        Eigen::VectorXd v(p);
        for (int k = 0; k != p; ++k)
            v[k] = sd * normal(engine);
        acc << v;
    }
    return acc.finalize();
}

} // namespace

TEST(t2_test, fewer_observations_than_components_one_sample)
{
    const Eigen::VectorXd zero = Eigen::VectorXd::Zero(40);
    const Eigen::MatrixXd l = mixing(40, 1.0);
    const int rejected = rejected_of_1000([&](std::mt19937 &engine) {
        return alps::alea::test_mean(sample(engine, zero, l, 20), zero);
    });
    EXPECT_GE(rejected, 20);
    EXPECT_LE(rejected, 100);
}

TEST(t2_test, fewer_observations_than_components_two_sample)
{
    const Eigen::VectorXd zero = Eigen::VectorXd::Zero(40);
    const Eigen::MatrixXd wide = mixing(40, 3.0), narrow = mixing(40, 1.0);
    const int rejected = rejected_of_1000([&](std::mt19937 &engine) {
        return alps::alea::test_mean(sample(engine, zero, wide, 20),
                                     sample(engine, zero, narrow, 30));
    });
    EXPECT_GE(rejected, 20);
    EXPECT_LE(rejected, 100);
}

TEST(t2_test, fewer_observations_than_components_independent)
{
    const Eigen::VectorXd zero = Eigen::VectorXd::Zero(50);
    const int rejected = rejected_of_1000([&](std::mt19937 &engine) {
        return alps::alea::test_mean(independent_sample(engine, 50, 1.0, 20), zero);
    });
    EXPECT_GE(rejected, 20);
    EXPECT_LE(rejected, 100);
}

// The fallback detects a shift, keeps NaN for NaN data, and needs at least 15
// observations (14 degrees of freedom): with fewer it is NaN, as before.
TEST(t2_test, fewer_observations_than_components_limits)
{
    std::mt19937 engine(29);
    const Eigen::VectorXd zero = Eigen::VectorXd::Zero(40);
    const Eigen::MatrixXd l = mixing(40, 1.0);
    Eigen::VectorXd shifted = zero;
    shifted.head(10).array() += 1.5;   // 10 components shifted by ~1.5 sd
    cov_result<double> res = sample(engine, shifted, l, 20);
    alps::alea::t2_result r = alps::alea::test_mean(res, zero);
    EXPECT_TRUE(r.normal());
    EXPECT_LT(r.pvalue(), 1e-6);

    alps::alea::t2_result few = alps::alea::test_mean(sample(engine, zero, l, 10), zero);
    EXPECT_TRUE(std::isnan(few.pvalue()));

    Eigen::VectorXd nan_expected = zero;
    nan_expected[3] = NAN;
    EXPECT_TRUE(std::isnan(alps::alea::test_mean(res, nan_expected).pvalue()));
}
