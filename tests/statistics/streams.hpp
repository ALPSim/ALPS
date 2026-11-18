// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Deterministic input streams with known statistics, shared by the tests in
// this directory.
#ifndef ALPS_TESTS_STATISTICS_STREAMS_HPP
#define ALPS_TESTS_STATISTICS_STREAMS_HPP

#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

namespace statistics_test {

// Standard normal deviates by Box-Muller on mt19937_64's raw output, which the
// C++ standard fixes; std::normal_distribution differs between standard
// libraries.
class normal {
  public:
    explicit normal(std::uint64_t seed) : engine_(seed) {}
    double operator()() {
        const double pi = 3.14159265358979323846;
        const double u1 = (engine_() >> 11) * 0x1.0p-53;
        const double u2 = (engine_() >> 11) * 0x1.0p-53;
        return std::sqrt(-2 * std::log1p(-u1)) * std::cos(2 * pi * u2);
    }
    // Uniform deviate in [0, 1).
    double uniform() { return (engine_() >> 11) * 0x1.0p-53; }
  private:
    std::mt19937_64 engine_;
};

// The AR(1) process x_t = phi x_{t-1} + sqrt(1 - phi^2) xi_t, started from its
// stationary distribution: mean 0, variance 1, and integrated autocorrelation
// time tau = phi / (1 - phi) in the convention error^2 = (1 + 2 tau) var / N.
class ar1 {
  public:
    ar1(std::uint64_t seed, double phi) : xi_(seed), phi_(phi), x_(xi_()) {}
    double operator()() {
        const double value = x_;
        x_ = phi_ * x_ + std::sqrt(1 - phi_ * phi_) * xi_();
        return value;
    }
  private:
    normal xi_;
    double phi_, x_;
};

inline std::vector<double> ar1_series(std::uint64_t seed, double phi, std::size_t n) {
    ar1 source(seed, phi);
    std::vector<double> values(n);
    for (double & value : values)
        value = source();
    return values;
}

// Exact variance of the sum of n consecutive values of the AR(1) process.
inline double ar1_sum_variance(double phi, double n) {
    return n * (1 + phi) / (1 - phi) - 2 * phi * (1 - std::pow(phi, n)) / ((1 - phi) * (1 - phi));
}

// Exact standard deviation of the mean of n consecutive values.
inline double ar1_exact_error(double phi, double n) {
    return std::sqrt(ar1_sum_variance(phi, n)) / n;
}

inline double ar1_tau(double phi) { return phi / (1 - phi); }

// Expected ratio of the error estimated from bins of length n to the exact
// error, in the limit of many bins: bins shorter than a few autocorrelation
// times are still correlated, so the estimate is too small by this factor.
inline double ar1_bin_factor(double phi, double n) {
    return std::sqrt(1 - 2 * phi * (1 - std::pow(phi, n)) / (n * (1 - phi * phi)));
}

// Textbook standard error of the mean from contiguous bins of the given length;
// a remainder shorter than one bin is left out of the bins, not of the mean.
inline double binned_error(const std::vector<double> & values, std::size_t bin_length) {
    double sum = 0;
    for (double value : values)
        sum += value;
    const double mean = sum / values.size();
    const std::size_t bins = values.size() / bin_length;
    double squares = 0;
    for (std::size_t b = 0; b < bins; ++b) {
        double bin = 0;
        for (std::size_t i = 0; i < bin_length; ++i)
            bin += values[b * bin_length + i];
        bin /= bin_length;
        squares += (bin - mean) * (bin - mean);
    }
    return std::sqrt(squares / (double(bins) * (bins - 1)));
}

}

#endif
