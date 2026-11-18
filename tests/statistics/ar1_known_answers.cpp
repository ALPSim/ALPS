// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Error bars and autocorrelation times against the exact values of the AR(1)
// process. A binned error from bins of length n is expected at
// exact * ar1_bin_factor(phi, n): each library is held to the value its own bin
// length implies, and to a floor on how far below the exact error it may report.
#include "streams.hpp"

#include <alps/alea.h>
#include <alps/alea.hpp>
#include <alps/accumulators.hpp>

#include <gtest/gtest.h>

#include <cmath>

namespace {

using statistics_test::ar1;
using statistics_test::ar1_bin_factor;
using statistics_test::ar1_exact_error;
using statistics_test::ar1_tau;

const std::size_t n = std::size_t(1) << 18;
const int replicas = 8;

// A library may report an error below the exact one only by this factor.
const double calibration_floor = 0.9;

struct averages {
    double naive = 0, error = 0, tau = 0;
};

// Relative tolerance of an error averaged over the replicas, each estimated
// from the given number of bins: four standard deviations.
double error_tolerance(double bins) { return 4 / std::sqrt(2 * (bins - 1) * replicas); }

void expect_calibrated(const averages & avg, double phi, double bin_length) {
    const double exact = ar1_exact_error(phi, n);
    const double factor = ar1_bin_factor(phi, bin_length);
    const double tolerance = error_tolerance(n / bin_length);
    EXPECT_NEAR(factor, avg.error / exact, tolerance * factor) << "phi = " << phi;
    // The naive error ignores the correlation.
    EXPECT_NEAR(std::sqrt((1 - phi) / (1 + phi)), avg.naive / exact, 0.02) << "phi = " << phi;
    // tau follows from the ratio of the binned to the naive variance.
    const double tau = ((1 + 2 * ar1_tau(phi)) * factor * factor - 1) / 2;
    EXPECT_NEAR(tau, avg.tau, 2 * tolerance * (tau + 0.5) + 0.02) << "phi = " << phi;
}

averages alea(double phi) {
    averages avg;
    for (int r = 0; r < replicas; ++r) {
        alps::SimpleRealObservable simple("simple");
        alps::RealObservable binned("binned");
        ar1 source(1000003u * r + 17, phi);
        for (std::size_t t = 0; t < n; ++t) {
            const double x = source();
            simple << x;
            binned << x;
        }
        avg.naive += simple.error() / replicas;
        avg.error += binned.error() / replicas;
        avg.tau += binned.tau() / replicas;
    }
    return avg;
}

averages accumulators(double phi) {
    namespace acc = alps::accumulators;
    averages avg;
    for (int r = 0; r < replicas; ++r) {
        acc::accumulator_set set;
        set << acc::NoBinningAccumulator<double>("simple") << acc::LogBinningAccumulator<double>("binned");
        ar1 source(1000003u * r + 17, phi);
        for (std::size_t t = 0; t < n; ++t) {
            const double x = source();
            set["simple"] << x;
            set["binned"] << x;
        }
        const acc::result_set results(set);
        avg.naive += results["simple"].error<double>() / replicas;
        avg.error += results["binned"].error<double>() / replicas;
        avg.tau += results["binned"].autocorrelation<double>() / replicas;
    }
    return avg;
}

averages aleav2(double phi) {
    averages avg;
    for (int r = 0; r < replicas; ++r) {
        alps::alea::var_acc<double> simple(1);
        alps::alea::autocorr_acc<double> binned(1);
        ar1 source(1000003u * r + 17, phi);
        for (std::size_t t = 0; t < n; ++t) {
            const double x = source();
            simple << x;
            binned << x;
        }
        const alps::alea::autocorr_result<double> result = binned.finalize();
        avg.naive += simple.finalize().stderror()[0] / replicas;
        avg.error += result.stderror()[0] / replicas;
        avg.tau += result.tau()[0] / replicas;
    }
    return avg;
}

class Ar1KnownAnswers : public ::testing::TestWithParam<double> {};

}

// ALEA and alps::accumulators report the error from 128 to 255 bins: bins of
// 2048 samples here.
TEST_P(Ar1KnownAnswers, Alea)
{
    const double phi = GetParam();
    const averages avg = alea(phi);
    expect_calibrated(avg, phi, 2048);
    EXPECT_GE(avg.error / ar1_exact_error(phi, n), calibration_floor);
}

TEST_P(Ar1KnownAnswers, Accumulators)
{
    const double phi = GetParam();
    const averages avg = accumulators(phi);
    expect_calibrated(avg, phi, 2048);
    EXPECT_GE(avg.error / ar1_exact_error(phi, n), calibration_floor);
}

// ALEA v2 reports the error from the highest level with at least 128 batches:
// batches of 2048 samples here.
TEST_P(Ar1KnownAnswers, AleaV2)
{
    const double phi = GetParam();
    const averages avg = aleav2(phi);
    expect_calibrated(avg, phi, 2048);
    EXPECT_GE(avg.error / ar1_exact_error(phi, n), calibration_floor);
}

INSTANTIATE_TEST_SUITE_P(Phi, Ar1KnownAnswers, ::testing::Values(0.0, 0.5, 0.9, 0.99));

// A run of about a hundred autocorrelation times cannot give a trustworthy
// error; a library must say so. 1000 samples at phi = 0.9: exact error 0.137.
TEST(ShortRun, AleaDoesNotClaimConvergence)
{
    alps::RealObservable binned("binned");
    ar1 source(5, 0.9);
    for (int t = 0; t < 1000; ++t)
        binned << source();
    EXPECT_NE(binned.converged_errors(), alps::CONVERGED);
}

TEST(ShortRun, AleaV2DoesNotClaimConvergence)
{
    alps::alea::autocorr_acc<double> binned(1);
    ar1 source(5, 0.9);
    for (int t = 0; t < 1000; ++t)
        binned << source();
    const alps::alea::autocorr_result<double> result = binned.finalize();
    EXPECT_NE(result.converged()[0], alps::alea::CONVERGED);
    // 250 batches of 4 samples: more than the naive error (0.031), still short
    // of the exact one, and the flag says so.
    EXPECT_GT(result.stderror()[0], 1.5 * result.level(0).stderror()[0]);
    EXPECT_GT(result.tau()[0], 1);
}

// Below 256 samples no level above the lowest has 128 batches. Nothing is
// known about the correlation, and tau must not read as zero.
TEST(ShortRun, AleaV2HasNoTauWithoutBatches)
{
    alps::alea::autocorr_acc<double> binned(1);
    ar1 source(5, 0.9);
    for (int t = 0; t < 200; ++t)
        binned << source();
    const alps::alea::autocorr_result<double> result = binned.finalize();
    EXPECT_FALSE(result.tau_available());
    EXPECT_TRUE(std::isnan(result.tau()[0]));
    EXPECT_NE(result.converged()[0], alps::alea::CONVERGED);
    // The minimum is a parameter: with 16 batches, level 3 qualifies.
    EXPECT_TRUE(result.tau_available(16));
    EXPECT_EQ(result.level(3).stderror()[0], result.stderror(16)[0]);
}

// The convergence flag separates a long run from a short one.
TEST(ConvergenceFlag, Alea)
{
    alps::RealObservable long_run("long"), short_run("short");
    ar1 fast(11, 0.5), slow(12, 0.99);
    for (std::size_t t = 0; t < n; ++t)
        long_run << fast();
    for (std::size_t t = 0; t < (std::size_t(1) << 14); ++t)
        short_run << slow();
    EXPECT_EQ(long_run.converged_errors(), alps::CONVERGED);
    EXPECT_NE(short_run.converged_errors(), alps::CONVERGED);
}

TEST(ConvergenceFlag, AleaV2)
{
    alps::alea::autocorr_acc<double> long_run(1), short_run(1);
    ar1 fast(11, 0.5), slow(12, 0.99);
    for (std::size_t t = 0; t < n; ++t)
        long_run << fast();
    for (std::size_t t = 0; t < (std::size_t(1) << 14); ++t)
        short_run << slow();
    EXPECT_EQ(long_run.finalize().converged()[0], alps::alea::CONVERGED);
    EXPECT_EQ(short_run.finalize().converged()[0], alps::alea::NOT_CONVERGED);
}
