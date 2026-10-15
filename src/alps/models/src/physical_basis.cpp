// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/src/operators/physical_basis.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports;
// HDF5 save/load removed.
#include <alps/models/physical_basis.hpp>
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace alps::models {
namespace {
void require(bool condition, std::string const& message) {
    if (!condition) throw std::invalid_argument(message);
}
} // namespace

void validate_basis(PhysicalBasis const& basis) {
    auto const& names = basis.conserved_quantum_numbers;
    require(std::set<std::string>(names.begin(), names.end()).size() == names.size(),
            "conserved quantum-number names must be distinct");
    require(basis.fermion_convention == FermionConvention::strict
            || basis.fermion_convention == FermionConvention::lattice,
            "fermion convention is unknown");
    for (auto const& name : names) require(!name.empty(), "quantum-number names must not be empty");
    for (auto const& site : basis.sites) {
        require(!site.labels.empty(), "physical basis must not be empty");
        require(site.statistics == Statistics::bosonic || site.statistics == Statistics::fermionic,
                "exchange statistics are unknown");
        require(site.parity.size() == site.labels.size(), "parity must match the physical basis dimension");
        for (double value : site.parity)
            require(value == 1. || value == -1., "parity values must be +1 or -1");
        for (auto const& [name, values] : site.quantum_numbers) {
            require(!name.empty(), "quantum-number names must not be empty");
            require(values.size() == site.labels.size(), "quantum-number labels must match the physical basis dimension");
            require(std::all_of(values.begin(), values.end(), [](double v) { return std::isfinite(v); }),
                    "quantum-number labels must be finite");
        }
        for (auto const& name : names)
            require((site.quantum_numbers.count(name) != 0), "physical basis lacks conserved quantum number " + name);
    }
}

void require_same_basis(PhysicalBasis const& actual, PhysicalBasis const& expected) {
    validate_basis(actual);
    validate_basis(expected);
    require(actual.sites.size() == expected.sites.size(), "physical basis lengths differ");
    require(actual.conserved_quantum_numbers == expected.conserved_quantum_numbers,
            "conserved quantum-number names or order differ");
    require(actual.fermion_convention == expected.fermion_convention,
            "fermion convention differs");
    for (std::size_t i = 0; i < actual.sites.size(); ++i) {
        auto const& a = actual.sites[i];
        auto const& e = expected.sites[i];
        auto const where = " at site " + std::to_string(i);
        require(a.labels.size() == e.labels.size(), "physical basis dimension differs" + where);
        require(a.labels == e.labels, "physical basis labels or order differ" + where);
        require(a.statistics == e.statistics, "exchange statistics differ" + where);
        require(a.parity == e.parity, "physical basis parity differs" + where);
        require(a.quantum_numbers == e.quantum_numbers, "physical basis quantum numbers differ" + where);
    }
}

} // namespace alps::models
