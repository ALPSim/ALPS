// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Shared helpers of the classic DMRG tests: exact-diagonalization references
// through ALPS::ed and DMRG options that keep enough states for the small
// systems under test to be solved exactly.
#pragma once

#include <alps/dmrg/classic.hpp>
#include <alps/ed/product_problem.hpp>
#include <alps/models/builtin.hpp>

#include <complex>
#include <map>
#include <string>
#include <vector>

namespace lat = alps::lattice;
namespace mdl = alps::models;
namespace cd = alps::dmrg::classic;

namespace classic_test {

using Sector = std::map<std::string, double>;

inline alps::ed::ProductEigenResult ed_states(mdl::Hamiltonian const& model, Sector const& sector, int count = 1) {
    alps::ed::ProductOptions options;
    options.sector = sector;
    alps::ed::EigenstateOptions request;
    request.count = count;
    request.method = alps::ed::EigenstateMethod::dense;
    return alps::ed::prepare_product(model, options).eigenstates(request);
}

inline std::vector<double> ed_energies(mdl::Hamiltonian const& model, Sector const& sector, int count = 1) {
    return ed_states(model, sector, count).eigenvalues;
}

// <psi_0|term|psi_0> for every term of `op`, as the classic solver reports it.
inline std::vector<double> ed_term_values(alps::ed::ProductEigenResult const& states,
                                          mdl::CompiledOperator const& op, std::size_t state = 0) {
    std::vector<double> values;
    for (auto const& term : op.terms) {
        mdl::CompiledOperator single{op.site_dims, {term}};
        values.push_back(alps::ed::expectation_values(states, single).at(state).real());
    }
    return values;
}

inline double ed_value(alps::ed::ProductEigenResult const& states, mdl::CompiledOperator const& op,
                       std::size_t state = 0) {
    return alps::ed::expectation_values(states, op).at(state).real();
}

// Enough kept states for an exact finite-system sweep on the small chains.
inline cd::Options exact_options(std::size_t states, Sector sector = {}, std::size_t count = 1) {
    cd::Options options;
    options.sweeps = cd::fixed_sweeps(4, states);
    options.warmup_states = states;
    options.num_eigenvalues = count;
    options.sector = std::move(sector);
    options.lanczos_tolerance = 1e-13;
    return options;
}

} // namespace classic_test
