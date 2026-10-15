// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/src/operators/basis.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
#include <alps/models/hamiltonian.hpp>
#include <alps/models/site.hpp>

#include <algorithm>
#include <stdexcept>

namespace alps::models {
namespace {
LocalBasis local_basis(SiteBasisPtr const& site) {
    if (!site || site->dim <= 0)
        throw std::invalid_argument("physical basis conversion requires positive site dimensions");
    if (site->statistics != Statistics::bosonic && site->statistics != Statistics::fermionic)
        throw std::invalid_argument("unknown model exchange statistics");
    auto labels = site->basis_labels;
    if (labels.empty()) {
        labels.reserve(static_cast<std::size_t>(site->dim));
        for (int state = 0; state < site->dim; ++state)
            labels.push_back(std::to_string(state));
    }
    if (labels.size() != static_cast<std::size_t>(site->dim))
        throw std::invalid_argument("physical basis labels must match the local dimension");
    return {std::move(labels), site->statistics, site->quantum_numbers, site->parity_operator};
}

void resolve_neutral_rows(LocalBasis& site, std::vector<std::string> const& names) {
    for (auto const& name : names)
        if (!(site.quantum_numbers.count(name) != 0))
            site.quantum_numbers.emplace(name, std::vector<double>(site.labels.size(), 0.));
}
} // namespace

PhysicalBasis physical_basis(std::vector<SiteBasisPtr> const& sites,
                          std::vector<std::string> const& conserved_qns,
                          FermionConvention convention) {
    if (convention != FermionConvention::strict && convention != FermionConvention::lattice)
        throw std::invalid_argument("unknown model fermion convention");
    PhysicalBasis result{{}, conserved_qns, convention};
    result.sites.reserve(sites.size());
    for (auto const& site : sites) result.sites.push_back(local_basis(site));
    for (auto const& name : conserved_qns) {
        if (std::none_of(result.sites.begin(), result.sites.end(), [&](auto const& site) {
                return (site.quantum_numbers.count(name) != 0);
            }))
            throw std::invalid_argument("no site defines declared conserved quantum number " + name);
    }
    for (auto& site : result.sites) resolve_neutral_rows(site, conserved_qns);
    validate_basis(result);
    return result;
}

PhysicalBasis physical_basis(Hamiltonian const& model) {
    std::vector<SiteBasisPtr> sites;
    sites.reserve(static_cast<std::size_t>(model.num_sites()));
    for (int site = 0; site < model.num_sites(); ++site)
        sites.push_back(model.site_for(site));
    return physical_basis(sites, model.conserved_qns(), model.fermion_convention());
}

void require_compatible_sites(PhysicalBasis const& reference,
                             std::vector<std::pair<std::size_t, SiteBasisPtr>> const& sites) {
    validate_basis(reference);
    for (auto const& [index, site] : sites) {
        if (index >= reference.sites.size())
            throw std::invalid_argument("operator site is outside the physical basis");
        auto actual = local_basis(site);
        resolve_neutral_rows(actual, reference.conserved_quantum_numbers);
        if (actual != reference.sites[index])
            throw std::invalid_argument("operator uses a different physical basis at site " + std::to_string(index));
    }
}

} // namespace alps::models
