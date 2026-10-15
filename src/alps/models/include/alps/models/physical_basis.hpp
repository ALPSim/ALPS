// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/include/alps/model/physical_basis.hpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports;
// HDF5 save/load removed.
#pragma once

#include <alps/models_export.h>
#include <map>
#include <string>
#include <vector>

namespace alps::models {

enum class Statistics { bosonic, fermionic };
enum class FermionConvention { strict, lattice };

// Entries in every array refer to the same ordered physical basis. Quantum
// numbers include all declared labels, including those not conserved by a run.
struct LocalBasis {
    std::vector<std::string> labels;
    Statistics statistics = Statistics::bosonic;
    std::map<std::string, std::vector<double>> quantum_numbers;
    std::vector<double> parity;
    bool operator==(LocalBasis const& o) const {
        return labels == o.labels && statistics == o.statistics
            && quantum_numbers == o.quantum_numbers && parity == o.parity;
    }
    bool operator!=(LocalBasis const& o) const { return !(*this == o); }
};

// Ordered physical local spaces of a model, independent of any Hamiltonian,
// operator registry, product-state enumeration or solver representation.
struct PhysicalBasis {
    std::vector<LocalBasis> sites;
    // This order defines the components of solver charge vectors.
    std::vector<std::string> conserved_quantum_numbers;
    FermionConvention fermion_convention = FermionConvention::strict;
    bool operator==(PhysicalBasis const& o) const {
        return sites == o.sites && conserved_quantum_numbers == o.conserved_quantum_numbers
            && fermion_convention == o.fermion_convention;
    }
    bool operator!=(PhysicalBasis const& o) const { return !(*this == o); }
};

// The empty basis represents the scalar vacuum. Every occupied site requires
// explicit labels, finite quantum numbers and one +/-1 parity per state.
ALPS_MODELS_DECL void validate_basis(PhysicalBasis const& basis);
ALPS_MODELS_DECL void require_same_basis(PhysicalBasis const& actual, PhysicalBasis const& expected);

} // namespace alps::models
