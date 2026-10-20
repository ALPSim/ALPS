/*
 * Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
 */

// Error bars against exact values. The AR(1) process
//     x_t = phi x_{t-1} + sqrt(1 - phi^2) xi_t,   xi_t ~ N(0, 1),
// started from its stationary distribution, has unit variance, integrated
// autocorrelation time tau = phi / (1 - phi) in the accumulators' convention
// (error^2 = (1 + 2 tau) var / N), and a known variance of the sample mean.

#include "alps/accumulators.hpp"
#include "gtest/gtest.h"

#include <cmath>
#include <cstdint>
#include <random>

namespace acc = alps::accumulators;

namespace {

const double phi = 0.9;
const long samples = 1L << 18;
const int replicas = 8;

// Box-Muller on mt19937_64's raw output, which the C++ standard fixes;
// std::normal_distribution differs between standard libraries.
class normal {
  public:
    explicit normal(std::uint64_t seed) : engine_(seed) {}
    double operator()() {
        const double pi = 3.14159265358979323846;
        double u1 = (engine_() >> 11) * 0x1.0p-53;
        double u2 = (engine_() >> 11) * 0x1.0p-53;
        return std::sqrt(-2 * std::log1p(-u1)) * std::cos(2 * pi * u2);
    }
  private:
    std::mt19937_64 engine_;
};

// Exact standard deviation of the mean of n samples.
double exact_error(long n) {
    double var = ((1 + phi) / (1 - phi) - 2 * phi * (1 - std::pow(phi, n)) / (n * (1 - phi) * (1 - phi))) / n;
    return std::sqrt(var);
}

struct averages {
    double no_binning_error = 0, log_binning_error = 0, full_binning_error = 0;
    double log_binning_tau = 0, full_binning_tau = 0;
};

averages run() {
    averages avg;
    for (int r = 0; r < replicas; ++r) {
        acc::accumulator_set measurements;
        measurements << acc::NoBinningAccumulator<double>("no_binning")
                     << acc::LogBinningAccumulator<double>("log_binning")
                     << acc::FullBinningAccumulator<double>("full_binning");
        normal xi(1000003u * r + 17);
        double x = xi();
        for (long t = 0; t < samples; ++t) {
            measurements["no_binning"] << x;
            measurements["log_binning"] << x;
            measurements["full_binning"] << x;
            x = phi * x + std::sqrt(1 - phi * phi) * xi();
        }
        acc::result_set results(measurements);
        avg.no_binning_error += results["no_binning"].error<double>() / replicas;
        avg.log_binning_error += results["log_binning"].error<double>() / replicas;
        avg.full_binning_error += results["full_binning"].error<double>() / replicas;
        avg.log_binning_tau += results["log_binning"].autocorrelation<double>() / replicas;
        avg.full_binning_tau += results["full_binning"].autocorrelation<double>() / replicas;
    }
    return avg;
}

}

TEST(AccumulatorCalibration, AR1) {
    const averages avg = run();
    const double error = exact_error(samples);
    const double tau = phi / (1 - phi);

    // Binning recovers the exact error. Each replica's estimate scatters by
    // about 9%, from the 256 bins it uses, so the average of 8 by about 3%.
    EXPECT_NEAR(1.0, avg.log_binning_error / error, 0.12);
    EXPECT_NEAR(1.0, avg.full_binning_error / error, 0.12);
    EXPECT_NEAR(tau, avg.log_binning_tau, 0.25 * tau);
    EXPECT_NEAR(tau, avg.full_binning_tau, 0.25 * tau);

    // No binning ignores the correlation: its error is the naive
    // sqrt(var / N), smaller by sqrt((1 - phi) / (1 + phi)).
    EXPECT_NEAR(std::sqrt((1 - phi) / (1 + phi)), avg.no_binning_error / error, 0.02);
}
