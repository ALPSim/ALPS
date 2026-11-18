// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// The observable layer on ALEA v2, <alps/alea/observables.hpp>, held to the
// same answers as the libraries it is meant to replace.
#include "streams.hpp"

#include <alps/alea.hpp>
#include <alps/alea/observables.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {

using alps::alea::histogram;
using alps::alea::observable;
using alps::alea::observable_set;
using alps::alea::observable_result;
using alps::alea::result_set;
using statistics_test::ar1_exact_error;
using statistics_test::ar1_series;
using statistics_test::binned_error;

const std::size_t n = std::size_t(1) << 16;

// The textbook jackknife of f over equally long bins of several series.
struct jackknife {
    double mean, error;
};

template <typename F> jackknife hand_jackknife(const std::vector<std::vector<double> > & series, std::size_t bins, F f) {
    const std::size_t length = series[0].size() / bins;
    std::vector<double> total(series.size(), 0.);
    std::vector<std::vector<double> > bin(series.size(), std::vector<double>(bins, 0.));
    for (std::size_t s = 0; s < series.size(); ++s)
        for (std::size_t t = 0; t < bins * length; ++t) {
            bin[s][t / length] += series[s][t];
            total[s] += series[s][t];
        }
    std::vector<double> mean(series.size()), left(series.size());
    for (std::size_t s = 0; s < series.size(); ++s)
        mean[s] = total[s] / (bins * length);
    const double full = f(mean);
    double sum = 0, sum2 = 0;
    for (std::size_t b = 0; b < bins; ++b) {
        for (std::size_t s = 0; s < series.size(); ++s)
            left[s] = (total[s] - bin[s][b]) / ((bins - 1) * length);
        // The pseudo-value of bin b.
        const double value = bins * full - (bins - 1) * f(left);
        sum += value;
        sum2 += value * value;
    }
    const double average = sum / bins;
    return {average, std::sqrt((sum2 / bins - average * average) / (bins - 1))};
}

}

// One stream: the layer reports what the engine's binning analysis reports,
// and its bins are the hand-computed ones.
TEST(Observables, ScalarMatchesEngineAndHandBins)
{
    const std::vector<double> values = ar1_series(20261116, 0.9, n);
    observable_set set;
    set << observable("x") << observable("plain").without_bins();
    alps::alea::autocorr_acc<double> engine(1);
    for (double value : values) {
        set["x"] << value;
        set["plain"] << value;
        engine << value;
    }
    EXPECT_EQ(set["x"].count(), n);
    const result_set results = set.results();
    const observable_result x = results["x"], plain = results["plain"];
    const alps::alea::autocorr_result<double> reference = engine.finalize();

    EXPECT_EQ(x.count(), n);
    EXPECT_EQ(x.mean(), std::accumulate(values.begin(), values.end(), 0.) / n);
    EXPECT_EQ(x.error(), reference.stderror()[0]);
    EXPECT_NEAR(x.error(), binned_error(values, 512), 1e-10 * x.error());
    EXPECT_EQ(x.tau(), reference.tau()[0]);
    EXPECT_EQ(x.converged(), alps::alea::CONVERGED);
    EXPECT_NEAR(x.binning_errors()[4], binned_error(values, 16), 1e-10 * x.error());

    ASSERT_TRUE(x.has_bins());
    EXPECT_EQ(x.bin_size(), 512u);
    const std::vector<std::vector<double> > bins = x.bin_means();
    ASSERT_EQ(bins.size(), 128u);
    for (std::size_t i = 0; i < bins.size(); ++i)
        ASSERT_NEAR(bins[i][0], std::accumulate(values.begin() + 512 * i, values.begin() + 512 * (i + 1), 0.) / 512, 1e-13);

    // Without bins: the naive error, no autocorrelation time, no jackknife.
    EXPECT_EQ(plain.mean(), x.mean());
    EXPECT_NEAR(plain.error(), binned_error(values, 1), 1e-10 * plain.error());
    EXPECT_TRUE(std::isnan(plain.tau()));
    EXPECT_FALSE(plain.has_bins());
}

TEST(Observables, VectorLearnsItsSizeAndKeepsIt)
{
    observable_set set;
    set << observable("v").labels({"first", "second"});
    statistics_test::normal xi(3);
    std::vector<double> first, second;
    for (int t = 0; t < 4096; ++t) {
        first.push_back(xi());
        second.push_back(2 + 3 * xi());
        set["v"] << std::vector<double>({first.back(), second.back()});
    }
    EXPECT_THROW(set["v"] << 1.0, std::invalid_argument);

    const observable_result v = set.results()["v"];
    ASSERT_EQ(v.size(), 2u);
    EXPECT_EQ(v.labels()[1], "second");
    EXPECT_NEAR(v.mean_vector()[1], std::accumulate(second.begin(), second.end(), 0.) / 4096, 1e-13);
    EXPECT_NEAR(v.error_vector()[0], binned_error(first, 32), 1e-10);
    EXPECT_NEAR(v.error_vector()[1], binned_error(second, 32), 1e-10);
    EXPECT_THROW(v.mean(), std::logic_error);
}

// y = x + noise: the covariance of the means is ((1, 1), (1, 1.25)) / N, and
// its diagonal is the squared error.
TEST(Observables, CovarianceOfAVector)
{
    observable_set set;
    set << observable("v");
    statistics_test::normal xi(17);
    for (std::size_t t = 0; t < n; ++t) {
        const double x = xi();
        set["v"] << std::vector<double>({x, x + 0.5 * xi()});
    }
    const observable_result v = set.results()["v"];
    const std::vector<std::vector<double> > cov = v.covariance();
    EXPECT_NEAR(cov[0][0] * n, 1.0, 0.3);
    EXPECT_NEAR(cov[0][1] * n, 1.0, 0.3);
    EXPECT_NEAR(cov[1][1] * n, 1.25, 0.35);
    EXPECT_EQ(cov[0][1], cov[1][0]);
    EXPECT_NEAR(std::sqrt(cov[1][1]), v.error_vector()[1], 1e-12);
}

// U = <x^4> / <x^2>^2 = 3 with error sqrt(24 / N), as in nonlinear.cpp. The
// layer's answer is the textbook jackknife over its 128 bins, however the
// expression is written.
TEST(Observables, FunctionsAreOneJackknifeOverMeasuredBins)
{
    observable_set set;
    set << observable("A") << observable("B") << observable("plain").without_bins() << observable("rare");
    statistics_test::normal xi(424242);
    std::vector<std::vector<double> > series(2);
    for (std::size_t t = 0; t < n; ++t) {
        const double x = xi();
        series[0].push_back(x * x);
        series[1].push_back(x * x * x * x);
        set["A"] << series[0].back();
        set["B"] << series[1].back();
        set["plain"] << x;
        if (t % 3 == 0)
            set["rare"] << x;
    }
    const result_set results = set.results();
    const observable_result a = results["A"], b = results["B"];
    const jackknife expected =
        hand_jackknife(series, 128, [](const std::vector<double> & m) { return m[1] / (m[0] * m[0]); });

    const observable_result written = b / (a * a);
    const observable_result direct = alps::alea::transform(
        [](const std::vector<double> & m) { return std::vector<double>(1, m[1] / (m[0] * m[0])); }, {a, b});
    const observable_result stepwise = b / apply([](double m) { return m * m; }, a) * 1.0 + 0.0;
    for (const observable_result & u : {written, direct, stepwise}) {
        EXPECT_NEAR(u.mean(), expected.mean, 1e-11);
        EXPECT_NEAR(u.error(), expected.error, 1e-11);
    }
    EXPECT_NEAR(3.0, written.mean(), 4 * std::sqrt(24.0 / n));
    EXPECT_NEAR(1.0, written.error() / std::sqrt(24.0 / n), 0.25);

    // No bins, or bins that do not line up: an error, not a guess.
    EXPECT_THROW(a * results["plain"], std::invalid_argument);
    EXPECT_THROW(a * results["rare"], std::invalid_argument);
}

// <s O> / <s> = 1.5 with the error of signed_ratio.cpp.
TEST(Observables, SignedObservableIsARatioByJackknife)
{
    const double exact_error = std::sqrt((1.64 / 0.36 + 0.81 * 0.64 / 0.1296 - 2 * 0.9 * 0.64 / 0.216) / n);
    observable_set set;
    set << observable("Sign") << observable("O").signed_by("Sign");
    statistics_test::normal xi(8128);
    std::vector<std::vector<double> > series(2);
    for (std::size_t t = 0; t < n; ++t) {
        const double sign = xi.uniform() < 0.8 ? 1. : -1.;
        series[0].push_back(sign * (1 + 0.3 * sign + xi()));
        series[1].push_back(sign);
        set["Sign"] << sign;
        set["O"] << series[0].back();
    }
    const result_set results = set.results();
    const observable_result ratio = results["O"];
    const jackknife expected = hand_jackknife(series, 128, [](const std::vector<double> & m) { return m[0] / m[1]; });
    EXPECT_EQ(ratio.sign(), "Sign");
    EXPECT_NEAR(ratio.mean(), expected.mean, 1e-11);
    EXPECT_NEAR(ratio.error(), expected.error, 1e-11);
    EXPECT_NEAR(1.5, ratio.mean(), 4 * exact_error);
    EXPECT_NEAR(1.0, ratio.error() / exact_error, 0.25);

    // In a function, the signed observable is still numerator over sign.
    const jackknife squared =
        hand_jackknife(series, 128, [](const std::vector<double> & m) { return m[0] * m[0] / (m[1] * m[1]); });
    EXPECT_NEAR((ratio * ratio).error(), squared.error, 1e-11);

    observable_set orphan;
    orphan << observable("O").signed_by("Sign");
    orphan["O"] << 1.0;
    EXPECT_THROW(orphan.results(), std::runtime_error);
}

// Eight independent runs at phi = 0.9, merged: counts add up, no sample is
// lost, and the error is that of the pooled data (as ManyRuns in merge.cpp).
TEST(Observables, MergedRunsPoolTheirStatistics)
{
    const int runs = 8, replicas = 40;
    const std::size_t run_length = std::size_t(1) << 13;
    double error = 0, sum = 0;
    for (int r = 0; r < replicas; ++r) {
        observable_set pooled;
        sum = 0;
        for (int k = 0; k < runs; ++k) {
            observable_set run;
            run << observable("x") << histogram("h", -1, 1);
            for (double value : ar1_series(104729u * r + 31u * k + 7, 0.9, run_length)) {
                run["x"] << value;
                run["h"] << (value > 0 ? 1 : -1);
                sum += value;
            }
            pooled.merge(run.results());
        }
        const result_set results = pooled.results();
        ASSERT_EQ(results["x"].count(), runs * run_length);
        ASSERT_NEAR(results["x"].mean(), sum / (runs * run_length), 1e-13);
        ASSERT_EQ(results["h"].histogram()[0] + results["h"].histogram()[2], runs * run_length);
        ASSERT_EQ(results["h"].histogram()[1], 0u);
        error += results["x"].error() / replicas;
    }
    const double pooled_error = ar1_exact_error(0.9, run_length) / std::sqrt(double(runs));
    EXPECT_GE(error / pooled_error, 0.9);
    EXPECT_LE(error / pooled_error, 1.05);
}

TEST(Observables, DeclarationsAndReset)
{
    observable_set set;
    set << observable("x") << histogram("h", 0, 9);
    EXPECT_THROW(set << observable("x"), std::invalid_argument);
    EXPECT_THROW(set["y"], std::out_of_range);
    EXPECT_THROW(set["h"] << 10, std::out_of_range);
    EXPECT_THROW(set["h"] << 0.5, std::out_of_range);
    EXPECT_THROW(set.results()["x"].mean(), std::runtime_error);

    set["x"] << 1.0;
    set["h"] << 3;
    observable_set copy = set;
    set.reset();
    EXPECT_EQ(set["x"].count(), 0u);
    EXPECT_EQ(set["h"].count(), 0u);
    EXPECT_EQ(copy["x"].count(), 1u);
    // After a reset the size may be chosen again.
    set["x"] << std::vector<double>({1., 2.});
    EXPECT_EQ(set.results()["x"].size(), 2u);
    EXPECT_EQ(copy.results()["h"].histogram()[3], 1u);
}
