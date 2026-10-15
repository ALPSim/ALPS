// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, files
// src/solvers/diag/tests/diag_thermodynamics.cpp (a plain main, here as the
// Thermodynamics cases) and src/solvers/diag/tests/signed_thermodynamics.cpp.
// Adaptation: C++17 (std::span inputs replaced by std::vector), GoogleTest,
// <alps/ed/...> include names.
#include <alps/ed/projected_basis.hpp>
#include <alps/ed/projected_problem.hpp>
#include <alps/ed/spectrum.hpp>
#include <alps/ed/thermodynamics.hpp>
#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using alps::ed::thermal_surface;
using alps::ed::thermal_weights;
void close(double actual, double expected) {
    ASSERT_TRUE(std::isfinite(actual));
    EXPECT_NEAR(actual, expected, 2e-12);
}
void near(double actual, double expected) {
    ASSERT_TRUE(std::isfinite(actual));
    EXPECT_NEAR(actual, expected, 3e-12 * std::max(1., std::abs(expected)));
}
std::vector<double> const none{};
alps::ed::EigenstateOptions lanczos(int count) {
    alps::ed::EigenstateOptions o;
    o.count = count;
    o.method = alps::ed::EigenstateMethod::lanczos;
    return o;
}
} // namespace

TEST(Thermodynamics, TwoIndependentSpinsAtAReferenceField) {
    // Two independent spins 1/2 at a nonzero reference field. Their
    // partition function is [2 cosh(h/(2T))]^2, including degeneracy.
    std::vector<double> const energies{0.3, 0., 0., -0.3}, charges{-1., 0., 0., 1.};
    std::vector<double> const temperatures{0.2, 0.7, 2.}, fields{-1.1, 0., 0.4, 1.2};
    auto const surface = thermal_surface(energies, charges, temperatures, fields, 0.3);
    for (std::size_t t = 0; t < temperatures.size(); ++t)
        for (std::size_t f = 0; f < fields.size(); ++f) {
            double const beta = 1. / temperatures[t], h = fields[f], x = beta * h / 2.;
            double const q = std::tanh(x), chi = beta / (2. * std::pow(std::cosh(x), 2));
            double const free_energy = -2. * std::log(2. * std::cosh(x)) / beta;
            auto const& value = surface[t * fields.size() + f];
            close(value.energy, -h * q);
            close(value.free_energy, free_energy);
            close(value.entropy, beta * (-h * q - free_energy));
            close(value.specific_heat, beta * h * h * chi);
            close(value.charge, q);
            close(value.susceptibility, chi);
        }
    std::vector<double> const single{-0.75}, zero{0.}, one{1.};
    auto const value = thermal_surface(single, none, one, zero).front();
    close(value.energy, single[0]);
    close(value.free_energy, single[0]);
    close(value.entropy, 0.);
    close(value.specific_heat, 0.);
    close(value.charge, 0.);
    close(value.susceptibility, 0.);
}

TEST(Thermodynamics, ConstantOffsetsDoNotAlterFluctuations) {
    // Adding constants to energy or conserved charge cannot alter their
    // fluctuations. Raw second moments lose this property by cancellation.
    std::vector<double> const zero{0.};
    for (double offset : {0., 1e12, -1e12}) {
        std::vector<double> const levels{offset, offset + 1.}, q{offset, offset + 1.};
        for (double temperature : {0.125, 1., 10.}) {
            std::vector<double> const grid{temperature};
            auto const shifted = thermal_surface(levels, q, grid, zero).front();
            double const beta = 1. / temperature, p = 1. / (1. + std::exp(beta));
            close(shifted.energy, offset + p);
            close(shifted.charge, offset + p);
            close(shifted.specific_heat, beta * beta * p * (1. - p));
            close(shifted.susceptibility, beta * p * (1. - p));
            close(shifted.entropy, std::log1p(std::exp(-beta)) + beta * p);
        }
    }
}

TEST(Thermodynamics, ExtremeLevelsStayFinite) {
    std::vector<double> const zero{0.};
    // A frozen state has exactly zero entropy and fluctuations even when
    // squaring its absolute energy or charge would overflow.
    std::vector<double> const cold{1e-100}, very_negative{-1e200}, very_large{1e200};
    auto const frozen = thermal_surface(very_negative, very_large, cold, zero).front();
    close(frozen.energy, very_negative[0]);
    close(frozen.charge, very_large[0]);
    close(frozen.entropy, 0.);
    close(frozen.specific_heat, 0.);
    close(frozen.susceptibility, 0.);

    // The dimensionless fluctuation remains finite when raw E^2 overflows.
    std::vector<double> const wide_levels{0., 1e200}, hot{1e200};
    auto const wide = thermal_surface(wide_levels, none, hot, zero).front();
    double const p = 1. / (1. + std::exp(1.));
    close(wide.specific_heat, p * (1. - p));
    close(wide.entropy, std::log1p(std::exp(-1.)) + p);

    // Unoccupied high-energy levels must not turn 0 * infinity into NaN.
    auto const frozen_pair = thermal_surface(wide_levels, none, cold, zero).front();
    close(frozen_pair.energy, 0.);
    close(frozen_pair.entropy, 0.);
    close(frozen_pair.specific_heat, 0.);
}

TEST(Thermodynamics, InvalidInputsAreRejected) {
    std::vector<double> const energies{0.3, 0., 0., -0.3}, single{-0.75}, zero{0.}, one{1.};
    std::vector<double> const inf{std::numeric_limits<double>::infinity()};
    std::vector<double> const nan{std::numeric_limits<double>::quiet_NaN()};
    std::vector<double> const negative{-1.}, huge{std::numeric_limits<double>::max()};
    EXPECT_THROW(thermal_surface(none, none, one, zero), std::invalid_argument);
    EXPECT_THROW(thermal_surface(energies, one, one, zero), std::invalid_argument);
    EXPECT_THROW(thermal_surface(single, none, none, zero), std::invalid_argument);
    EXPECT_THROW(thermal_surface(single, none, one, none), std::invalid_argument);
    EXPECT_THROW(thermal_surface(inf, none, one, zero), std::invalid_argument);
    EXPECT_THROW(thermal_surface(nan, none, one, zero), std::invalid_argument);
    EXPECT_THROW(thermal_surface(single, inf, one, zero), std::invalid_argument);
    EXPECT_THROW(thermal_surface(single, none, zero, zero), std::invalid_argument);
    EXPECT_THROW(thermal_surface(single, none, nan, zero), std::invalid_argument);
    EXPECT_THROW(thermal_surface(single, none, inf, zero), std::invalid_argument);
    EXPECT_THROW(thermal_surface(single, none, one, inf), std::invalid_argument);
    EXPECT_THROW(thermal_surface(single, none, one, zero, inf[0]), std::invalid_argument);
    EXPECT_THROW(thermal_surface(negative, huge, one, huge), std::invalid_argument);
}

TEST(SignedThermodynamics, PositiveTwoLevelReferenceRetainsResponseAndOffsets) {
    std::vector<double> const energies{7., 9.}, charges{-3., -2.};
    std::vector<double> const temperatures{.25, 1., 4.}, fields{-.5, .25, 2.5};
    constexpr double field0 = .25;
    auto const result = thermal_surface(energies, charges, temperatures, fields, field0);
    for (std::size_t t = 0; t < temperatures.size(); ++t)
        for (std::size_t f = 0; f < fields.size(); ++f) {
            auto const temperature = temperatures[t];
            auto const shift = fields[f] - field0;
            auto const base = energies[0] - charges[0] * shift;
            auto const gap = 2. - shift;
            auto const x = gap / temperature;
            auto const occupied = 1. / (1. + std::exp(x));
            auto const& value = result[t * fields.size() + f];
            near(value.energy, base + gap * occupied);
            near(value.free_energy, base - temperature * std::log1p(std::exp(-x)));
            near(value.entropy, std::log1p(std::exp(-x)) + x * occupied);
            near(value.specific_heat, x * x * occupied * (1. - occupied));
            near(value.charge, charges[0] + occupied);
            near(value.susceptibility, occupied * (1. - occupied) / temperature);
        }
}

TEST(SignedThermodynamics, NegativeTemperatureInvertsPopulationAndResponse) {
    std::vector<double> const energies{0., 2.}, charges{1., 3.};
    std::vector<double> const temperatures{-2., .5, -1., 4.}, fields{0., .3, 2.};
    auto const result = thermal_surface(energies, charges, temperatures, fields);
    for (std::size_t t = 0; t < temperatures.size(); ++t)
        for (std::size_t f = 0; f < fields.size(); ++f) {
            auto const temperature = temperatures[t];
            auto const base = -fields[f];
            auto const gap = 2. - 2. * fields[f];
            auto const x = gap / temperature;
            auto const occupied = 1. / (1. + std::exp(x));
            auto const& value = result[t * fields.size() + f];
            near(value.energy, base + gap * occupied);
            near(value.free_energy, base - temperature * std::log1p(std::exp(-x)));
            near(value.entropy, std::log1p(std::exp(-x)) + x * occupied);
            near(value.specific_heat, x * x * occupied * (1. - occupied));
            near(value.charge, 1. + 2. * occupied);
            near(value.susceptibility, 4. * occupied * (1. - occupied) / temperature);
        }
}

TEST(SignedThermodynamics, SignedTemperaturesPreserveDegeneracyAndConstantEnergyShifts) {
    std::vector<double> const temperatures{-.25, .25, -8., 8.}, fields{0.};
    for (double offset : {0., 1e12, -1e12}) {
        std::vector<double> const levels{offset, offset, offset};
        auto const result = thermal_surface(levels, none, temperatures, fields);
        for (std::size_t t = 0; t < temperatures.size(); ++t) {
            near(result[t].energy, offset);
            near(result[t].free_energy, offset - temperatures[t] * std::log(3.));
            near(result[t].entropy, std::log(3.));
            near(result[t].specific_heat, 0.);
            near(result[t].partition, 3.);
        }
    }
}

TEST(SignedThermodynamics, FrozenInvertedEnsembleAvoidsOverflowAndUnoccupiedMoments) {
    std::vector<double> const levels{0., 1e200}, temperatures{-1e-100}, fields{0.};
    auto const value = thermal_surface(levels, none, temperatures, fields).front();
    near(value.energy, levels[1]);
    near(value.free_energy, levels[1]);
    near(value.entropy, 0.);
    near(value.specific_heat, 0.);
    near(value.partition, 1.);
}

TEST(SignedThermodynamics, ZeroAndNonfiniteTemperaturesRemainInvalid) {
    std::vector<double> const levels{0., 1.}, fields{0.};
    for (double temperature : {0., -0., std::numeric_limits<double>::infinity(),
                               -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::denorm_min()}) {
        std::vector<double> const grid{temperature};
        EXPECT_THROW(thermal_surface(levels, none, grid, fields), std::invalid_argument);
    }
}

TEST(SignedThermodynamics, PartitionReportsItsReferenceAndAbsoluteLogarithm) {
    std::vector<double> const fields{0.}, temperatures{-.5, 2.};
    for (double offset : {-1e12, 0., 1e12}) {
        std::vector<double> const levels{offset, offset + 2., offset + 2.};
        auto const result = thermal_surface(levels, none, temperatures, fields);
        for (std::size_t t = 0; t < temperatures.size(); ++t) {
            auto const temperature = temperatures[t];
            auto const reference = temperature > 0. ? levels[0] : levels[1];
            auto const expected = temperature > 0. ? 1. + 2. * std::exp(-2. / temperature) : 2. + std::exp(2. / temperature);
            EXPECT_DOUBLE_EQ(result[t].energy_reference, reference);
            near(result[t].partition, expected);
            near(result[t].log_partition, std::log(expected) - reference / temperature);
        }
    }
}

TEST(SignedThermodynamics, UnrepresentableLogPartitionDoesNotCorruptFiniteObservables) {
    std::vector<double> const fields{0.}, temperatures{-1e-300, 1e-300};
    for (double offset : {-1e300, 1e300}) {
        std::vector<double> const levels{offset, offset};
        auto const result = thermal_surface(levels, none, temperatures, fields);
        for (std::size_t t = 0; t < temperatures.size(); ++t) {
            auto const& value = result[t];
            near(value.energy, offset);
            near(value.free_energy, offset);
            near(value.entropy, std::log(2.));
            near(value.specific_heat, 0.);
            near(value.partition, 2.);
            ASSERT_TRUE(std::isinf(value.log_partition));
            EXPECT_EQ(std::signbit(value.log_partition), offset / temperatures[t] > 0.);
        }
    }
}

TEST(SignedThermodynamics, FiniteTemperatureScalesOverflowingEnergyDifferences) {
    constexpr double scale = 1e308;
    std::vector<double> const levels{-scale, scale}, fields{0.}, temperatures{-scale, scale};
    auto const result = thermal_surface(levels, none, temperatures, fields);
    auto const excited = 1. / (1. + std::exp(2.));
    for (std::size_t t = 0; t < temperatures.size(); ++t) {
        near(result[t].energy, -std::copysign(scale * std::tanh(1.), temperatures[t]));
        near(result[t].entropy, std::log1p(std::exp(-2.)) + 2. * excited);
        near(result[t].specific_heat, 4. * excited * (1. - excited));
        auto const weights = thermal_weights(levels, temperatures[t]);
        near(weights[t == 0 ? 0 : 1], excited);
        near(weights[0] + weights[1], 1.);
    }

    // The mean energy remains inside the finite energy bounds even when its
    // displacement from the chosen reference alone exceeds the double range.
    constexpr double wide = 1.7e308;
    std::vector<double> many(98, -wide);
    many.push_back(wide);
    std::vector<double> const inverted{-wide};
    auto const value = thermal_surface(many, none, inverted, fields).front();
    auto const low = 98. * std::exp(-2.) / (1. + 98. * std::exp(-2.));
    near(value.energy, wide * (1. - 2. * low));
}

TEST(SignedThermodynamics, ReusableWeightsMatchIndependentObservableAndRejectInvalidInputs) {
    std::vector<double> const levels{0., 2.};
    for (double temperature : {-2., .5}) {
        auto const weights = thermal_weights(levels, temperature);
        ASSERT_EQ(weights.size(), 2);
        auto const occupied = 1. / (1. + std::exp(2. / temperature));
        near(weights[1], occupied);
        near(-3. * weights[0] + 4. * weights[1], -3. + 7. * occupied);
    }
    std::vector<double> const nonfinite{std::numeric_limits<double>::infinity()};
    EXPECT_THROW(thermal_weights(none, 1.), std::invalid_argument);
    EXPECT_THROW(thermal_weights(nonfinite, 1.), std::invalid_argument);
    EXPECT_THROW(thermal_weights(levels, 0.), std::invalid_argument);
    EXPECT_THROW(thermal_weights(levels, std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
    std::vector<double> const equal{1e300, 1e300};
    auto const weights = thermal_weights(equal, -1e-300);
    EXPECT_DOUBLE_EQ(weights[0], .5);
    EXPECT_DOUBLE_EQ(weights[1], .5);
}

TEST(SignedThermodynamics, LargestFiniteTemperatureKeepsFiniteMeanDespiteReciprocalRounding) {
    double const scale = std::numeric_limits<double>::max();
    std::vector<double> const fields{0.};
    auto const minority = 98. * std::exp(-2.) / (1. + 98. * std::exp(-2.));
    for (double sign : {-1., 1.}) {
        std::vector<double> levels(98, sign * scale);
        levels.push_back(-sign * scale);
        std::vector<double> const temperatures{sign * scale};
        auto const value = thermal_surface(levels, none, temperatures, fields).front();
        near(value.energy, sign * scale * (2. * minority - 1.));

        // F itself remains finite although log(scaled_partition)/beta does
        // not. Cover the independent free-energy reconstruction fallback.
        std::vector<double> const equal(4, sign * scale);
        auto const degenerate = thermal_surface(equal, none, temperatures, fields).front();
        near(degenerate.free_energy, sign * scale * (1. - std::log(4.)));
    }
}

TEST(SignedThermodynamics, NativeSpectrumDimensionsDistinguishSelectionFromCompleteness) {
    std::vector<alps::ed::SpectrumLocalSpace> const spaces{{{1., 1.}, {}}};
    std::vector<std::int64_t> const offsets{0, 0, 1}, columns{1}, states{0, 1};
    std::vector<std::complex<double>> const values{2.};
    alps::models::CompiledOperator input;
    input.site_dims = {2};
    input.terms = {{1., {{0, {0., 0., 0., 2.}}}, {}}};
    for (int count : {1, 4}) {
        auto const csr = alps::ed::symmetry_spectrum(offsets, columns, values, 2, states, spaces, {}, count);
        auto const local = alps::ed::symmetry_spectrum_from_local(input, spaces, {}, {}, count);
        auto const basis = alps::ed::make_projected_basis(spaces, {});
        auto const projected = alps::ed::prepare_projected(input, basis).eigenstates(lanczos(count));
        ASSERT_EQ(csr.size(), 1);
        ASSERT_EQ(local.size(), 1);
        for (auto const* block : {&csr.front(), &local.front(), &projected.spectrum}) {
            EXPECT_EQ(block->block_dimension, 2);
            EXPECT_EQ(block->eigenvalues.size(), count == 1 ? 1 : 2);
        }
    }
}
