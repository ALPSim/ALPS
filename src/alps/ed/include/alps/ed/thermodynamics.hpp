// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/include/alps/diag/thermodynamics.hpp.
// Adaptation: C++17 (std::span replaced by std::vector const&), namespace alps::ed.
//
// Canonical-ensemble averages from a finite spectrum (k_B = 1), as the legacy
// fulldiag evaluation computes them from stored eigenvalues.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace alps::ed {

struct ThermalAverages {
    double energy, free_energy, entropy, specific_heat, charge, susceptibility;
    // Reference-scaled sum, not the absolute partition function:
    // sum(exp(-beta * (energy - energy_reference))).
    double partition;
    double energy_reference;
    // log(partition) - beta * energy_reference. This may be signed infinity
    // when the absolute log partition is outside the representable double range;
    // the other thermodynamic quantities do not depend on this representation.
    double log_partition;
};

namespace detail {
template<bool FiniteRange>
ThermalAverages accumulate_thermal_levels(std::vector<double> const&, std::vector<double> const&,
                                          double, double);
inline double inverse_temperature(double temperature) {
    if (!std::isfinite(temperature) || temperature == 0. || !std::isfinite(1. / temperature))
        throw std::invalid_argument("temperatures must be finite and nonzero with finite inverse");
    return 1. / temperature;
}

inline void validate_energies(std::vector<double> const& energies) {
    if (energies.empty()) throw std::invalid_argument("no eigenvalues: empty Hilbert space");
    for (double energy : energies)
        if (!std::isfinite(energy)) throw std::invalid_argument("eigenvalues must be finite");
}

inline double scaled_energy(double energy, double reference, double beta) {
    double difference = energy - reference;
    // Opposite finite extremes can overflow their difference even when its
    // dimensionless value is finite. Keep the usual arithmetic order otherwise.
    return std::isfinite(difference) ? beta * difference : beta * energy - beta * reference;
}
} // namespace detail

// The energy reference keeps the Boltzmann weights bounded. add() returns
// the same weight for optional state-resolved observable accumulation.
// Callers supply finite nonzero beta and finite levels/charges, with reference
// the lowest energy for beta > 0 and highest energy for beta < 0.
class ThermalAccumulator {
public:
    ThermalAccumulator(double beta, double reference) : beta_(beta), reference_(reference) {}

    double add(double energy, double charge) {
        return add_scaled(detail::scaled_energy(energy, reference_, beta_), charge);
    }

private:
    template<bool FiniteRange> friend ThermalAverages detail::accumulate_thermal_levels(
        std::vector<double> const&, std::vector<double> const&, double, double);

    double add_scaled(double scaled_energy, double charge) {
        double weight = std::exp(-scaled_energy);
        // Exactly unoccupied levels contribute nothing, even if their energy
        // separation overflows when scaled by the inverse temperature.
        if (weight == 0.) return weight;
        if (partition_ == 0.) charge_reference_ = charge;
        double relative_charge = charge - charge_reference_;
        partition_ += weight;
        energy_ += scaled_energy * weight;
        energy2_ += scaled_energy * scaled_energy * weight;
        charge_ += relative_charge * weight;
        charge2_ += relative_charge * relative_charge * weight;
        return weight;
    }

public:
    [[nodiscard]] ThermalAverages result() const {
        double mean_energy = energy_ / partition_;
        double log_partition = std::log(partition_);
        double mean_charge = charge_ / partition_;
        double relative_energy = mean_energy / beta_;
        // Halving before the fused reconstruction also handles finite temperatures
        // whose rounded inverse has an unrepresentable reciprocal (e.g. DBL_MAX).
        double energy = std::isfinite(relative_energy) ? reference_ + relative_energy
            : 2. * std::fma(mean_energy, .5 / beta_, .5 * reference_);
        double relative_free_energy = log_partition / beta_;
        double free_energy = std::isfinite(relative_free_energy) ? reference_ - relative_free_energy
            : 2. * std::fma(-log_partition, .5 / beta_, .5 * reference_);
        return {energy, free_energy,
                log_partition + mean_energy, energy2_ / partition_ - mean_energy * mean_energy,
                charge_reference_ + mean_charge,
                (charge2_ / partition_ - mean_charge * mean_charge) * beta_, partition_,
                reference_, std::fma(-beta_, reference_, log_partition)};
    }

private:
    double beta_, reference_;
    double charge_reference_ = 0.;
    // Energy moments use beta * (E - reference), so specific heat and
    // entropy never subtract large absolute energies or square them.
    double partition_ = 0., energy_ = 0., energy2_ = 0., charge_ = 0., charge2_ = 0.;
};

namespace detail {
// Keep ordinary and overflow-protected loops separate: carrying both sets of
// temporaries through the hot loop spills moment accumulators on arm64.
// FiniteRange is selected only after checking the complete energy span.
template<bool FiniteRange>
inline ThermalAverages accumulate_thermal_levels(std::vector<double> const& energies,
    std::vector<double> const& charges, double beta, double reference) {
    ThermalAccumulator accumulator(beta, reference);
    for (std::size_t i = 0; i < energies.size(); ++i) {
        if constexpr (FiniteRange)
            accumulator.add_scaled(beta * (energies[i] - reference), charges.empty() ? 0. : charges[i]);
        else
            accumulator.add(energies[i], charges.empty() ? 0. : charges[i]);
    }
    return accumulator.result();
}
} // namespace detail

// Normalized probabilities for a finite spectrum at a finite, nonzero signed
// temperature (k_B = 1). Per-state observables can be computed once and then
// contracted with these probabilities for each temperature without solving
// the Hamiltonian again. The supplied levels define the selected ensemble.
inline std::vector<double> thermal_weights(std::vector<double> const& energies, double temperature) {
    detail::validate_energies(energies);
    double beta = detail::inverse_temperature(temperature);
    double reference = beta > 0. ? *std::min_element(energies.begin(), energies.end())
                                 : *std::max_element(energies.begin(), energies.end());
    std::vector<double> weights;
    weights.reserve(energies.size());
    double partition = 0.;
    for (double energy : energies) {
        double weight = std::exp(-detail::scaled_energy(energy, reference, beta));
        weights.push_back(weight);
        partition += weight;
    }
    for (double& weight : weights) weight /= partition;
    return weights;
}

// Return a temperature-major surface (index t * fields.size() + f). charges is
// either empty (zero charge) or carries one conserved charge per eigenvalue.
// The field shifts each level by -charge * (field - field0), without another
// diagonalization. Temperatures are finite and nonzero; negative temperatures
// describe the population-inverted ensemble of the supplied finite, bounded
// spectrum.
inline std::vector<ThermalAverages> thermal_surface(
    std::vector<double> const& energies, std::vector<double> const& charges,
    std::vector<double> const& temperatures, std::vector<double> const& fields,
    double field0 = 0.) {
    detail::validate_energies(energies);
    if (!charges.empty() && charges.size() != energies.size())
        throw std::invalid_argument("charges must have one entry per eigenvalue");
    if (temperatures.empty() || fields.empty())
        throw std::invalid_argument("temperature and field grids must not be empty");
    for (double value : charges)
        if (!std::isfinite(value)) throw std::invalid_argument("charges must be finite");
    for (double value : temperatures)
        detail::inverse_temperature(value);
    for (double value : fields)
        if (!std::isfinite(value)) throw std::invalid_argument("fields must be finite");
    if (!std::isfinite(field0)) throw std::invalid_argument("reference field must be finite");

    std::vector<ThermalAverages> surface;
    if (temperatures.size() > surface.max_size() / fields.size())
        throw std::length_error("thermodynamic surface is too large");
    surface.resize(temperatures.size() * fields.size());
    std::vector<double> shifted(energies.size());
    for (std::size_t f = 0; f < fields.size(); ++f) {
        double minimum = std::numeric_limits<double>::infinity();
        double maximum = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < energies.size(); ++i) {
            double charge = charges.empty() ? 0. : charges[i];
            shifted[i] = energies[i] - charge * (fields[f] - field0);
            if (!std::isfinite(shifted[i]))
                throw std::invalid_argument("field-shifted eigenvalues must be finite");
            minimum = std::min(minimum, shifted[i]);
            maximum = std::max(maximum, shifted[i]);
        }
        for (std::size_t t = 0; t < temperatures.size(); ++t) {
            double beta = 1. / temperatures[t];
            double reference = beta > 0. ? minimum : maximum;
            surface[t * fields.size() + f] = std::isfinite(maximum - minimum)
                ? detail::accumulate_thermal_levels<true>(shifted, charges, beta, reference)
                : detail::accumulate_thermal_levels<false>(shifted, charges, beta, reference);
        }
    }
    return surface;
}

} // namespace alps::ed
