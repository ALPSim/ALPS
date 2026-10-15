// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/tests/models/heterogeneous_physical_basis.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports;
// MPO and HDF5 cases removed.
#include <alps/models/builtin.hpp>
#include <alps/models/physical_basis.hpp>
#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>

namespace {
namespace mdl = alps::models;

alps::lattice::Lattice kondo_cell() {
    alps::lattice::UnitCell cell(1);
    cell.add_vertex({0.}, 0);
    cell.add_vertex({.5}, 1);
    cell.add_edge(0, 1, 1);
    return cell.tile({1});
}

mdl::Hamiltonian explicit_kondo() {
    auto fermion = mdl::spinful_fermion();
    mdl::SiteBasis spin = *mdl::spin();
    spin.quantum_numbers["N"] = {0., 0.};
    auto neutral = mdl::finalize_site(std::move(spin));
    mdl::Hamiltonian model(kondo_cell(), {{0, fermion}, {1, neutral}},
                           mdl::FermionConvention::strict, {"N", "Sz"});
    auto f = [&](char const* name) { return mdl::TermSum::located(fermion, name, 0); };
    auto s = [&](char const* name) { return mdl::TermSum::located(neutral, name, 1); };
    model.add(f("Sz") * s("Sz") + .5 * (f("Sp") * s("Sm") + f("Sm") * s("Sp")));
    return model;
}

mdl::PhysicalBasis explicit_basis() {
    return {{{{"empty", "up", "down", "double"}, mdl::Statistics::fermionic,
              {{"N", {0., 1., 1., 2.}}, {"Sz", {0., .5, -.5, 0.}}}, {1., -1., -1., 1.}},
             {{"up", "down"}, mdl::Statistics::bosonic,
              {{"N", {0., 0.}}, {"Sz", {.5, -.5}}}, {1., 1.}}},
            {"N", "Sz"}, mdl::FermionConvention::strict};
}

} // namespace

TEST(HeterogeneousPhysicalBasis, NativeConversionPreservesAllDeclaredMeaningAndNeedsNoOperators) {
    auto fermion = std::make_shared<mdl::SiteBasis>(*mdl::spinful_fermion());
    auto spin = std::make_shared<mdl::SiteBasis>(*mdl::spin());
    fermion->basis_labels = {"empty", "up", "down", "double"};
    spin->basis_labels = {"up", "down"};
    spin->quantum_numbers["orbital"] = {2., 3.};
    // Basis conversion never needs identity matrices or an operator registry.
    fermion->operators.clear();
    spin->operators.clear();
    auto const converted = mdl::physical_basis({fermion, spin}, {"Sz", "N"}, mdl::FermionConvention::lattice);
    EXPECT_EQ(converted.conserved_quantum_numbers, (std::vector<std::string>{"Sz", "N"}));
    EXPECT_EQ(converted.fermion_convention, mdl::FermionConvention::lattice);
    EXPECT_EQ(converted.sites[1].labels, spin->basis_labels);
    EXPECT_EQ(converted.sites[1].quantum_numbers.at("orbital"), (std::vector<double>{2., 3.}));
    EXPECT_FALSE(converted.sites[0].quantum_numbers.count("orbital") != 0);
    EXPECT_FALSE(spin->quantum_numbers.count("N") != 0);
    spin->quantum_numbers["Sz"][0] = 19.;
    fermion->parity_operator[1] = 1.;
    EXPECT_EQ(converted.sites[1].quantum_numbers.at("Sz"), (std::vector<double>{.5, -.5}));
    EXPECT_EQ(converted.sites[0].parity[1], -1.);
    EXPECT_EQ(mdl::physical_basis(std::vector<mdl::SiteBasisPtr>{}), mdl::PhysicalBasis{});
}

TEST(HeterogeneousPhysicalBasis, NativeConversionRejectsUnknownOrMalformedMeaning) {
    auto spin = mdl::spin();
    EXPECT_THROW(mdl::physical_basis({spin}, {"N"}), std::invalid_argument);
    EXPECT_THROW(mdl::physical_basis({spin}, {"Sz", "Sz"}), std::invalid_argument);
    EXPECT_THROW(mdl::physical_basis({nullptr}, {"Sz"}), std::invalid_argument);
    EXPECT_THROW(mdl::physical_basis({spin}, {}, static_cast<mdl::FermionConvention>(99)), std::invalid_argument);
    auto changed = std::make_shared<mdl::SiteBasis>(*spin);
    changed->quantum_numbers["N"] = {0.};
    EXPECT_THROW(mdl::physical_basis({changed}, {"N"}), std::invalid_argument);
    changed->quantum_numbers["N"] = {0., std::numeric_limits<double>::quiet_NaN()};
    EXPECT_THROW(mdl::physical_basis({changed}, {"N"}), std::invalid_argument);
    changed->quantum_numbers.erase("N");
    changed->basis_labels = {"only one"};
    EXPECT_THROW(mdl::physical_basis({changed}), std::invalid_argument);
    changed->basis_labels.clear();
    changed->statistics = static_cast<mdl::Statistics>(99);
    EXPECT_THROW(mdl::physical_basis({changed}), std::invalid_argument);
}

TEST(HeterogeneousPhysicalBasis, OriginalOperatorSitesMatchTheirCompleteReferenceInOneBatch) {
    auto const model = mdl::kondo_lattice(kondo_cell(), 0., 1.);
    auto const reference = mdl::physical_basis(model);
    auto const original = reference;
    auto const spin = model.site_for(1);
    // This subset contains no model site declaring N. Its complete reference
    // supplies the charge context; a standalone full model still must declare N.
    EXPECT_NO_THROW(mdl::require_compatible_sites(reference, {{1, spin}}));
    EXPECT_NO_THROW(mdl::require_compatible_sites(reference,
        {{0, model.site_for(0)}, {1, spin}, {1, spin}}));
    auto explicit_spin = std::make_shared<mdl::SiteBasis>(*spin);
    explicit_spin->quantum_numbers["N"] = {0., 0.};
    EXPECT_NO_THROW(mdl::require_compatible_sites(reference, {{1, explicit_spin}}));
    EXPECT_THROW(mdl::physical_basis({spin}, {"N", "Sz"}), std::invalid_argument);
    EXPECT_FALSE(spin->quantum_numbers.count("N") != 0);
    EXPECT_EQ(reference, original);
}

TEST(HeterogeneousPhysicalBasis, ContextualSiteMatchingPreservesAllOtherPhysicalMetadata) {
    auto const model = mdl::kondo_lattice(kondo_cell(), 0., 1.);
    auto const reference = mdl::physical_basis(model);
    auto const spin = model.site_for(1);
    auto changed = std::make_shared<mdl::SiteBasis>(*spin);
    changed->basis_labels = {"1", "0"};
    EXPECT_THROW(mdl::require_compatible_sites(reference, {{1, changed}}), std::invalid_argument);
    *changed = *spin;
    changed->quantum_numbers["N"] = {1., 0.};
    EXPECT_THROW(mdl::require_compatible_sites(reference, {{1, changed}}), std::invalid_argument);
    *changed = *spin;
    changed->quantum_numbers["Sz"] = {1., -1.};
    EXPECT_THROW(mdl::require_compatible_sites(reference, {{1, changed}}), std::invalid_argument);
    *changed = *spin;
    changed->quantum_numbers["extra"] = {0., 0.};
    EXPECT_THROW(mdl::require_compatible_sites(reference, {{1, changed}}), std::invalid_argument);
    auto with_extra = reference;
    with_extra.sites[1].quantum_numbers["extra"] = {0., 0.};
    EXPECT_THROW(mdl::require_compatible_sites(with_extra, {{1, spin}}), std::invalid_argument);
    EXPECT_NO_THROW(mdl::require_compatible_sites(with_extra, {{1, changed}}));
    *changed = *spin;
    changed->statistics = mdl::Statistics::fermionic;
    EXPECT_THROW(mdl::require_compatible_sites(reference, {{1, changed}}), std::invalid_argument);
    *changed = *spin;
    changed->parity_operator[1] = -1.;
    EXPECT_THROW(mdl::require_compatible_sites(reference, {{1, changed}}), std::invalid_argument);
    EXPECT_THROW(mdl::require_compatible_sites(reference, {{2, spin}}), std::invalid_argument);
    EXPECT_THROW(mdl::require_compatible_sites(reference, {{1, nullptr}}), std::invalid_argument);
}

TEST(HeterogeneousPhysicalBasis, ContextualMatchingNeverRepairsAnIncompleteReference) {
    auto const model = mdl::kondo_lattice(kondo_cell(), 0., 1.);
    auto reference = mdl::physical_basis(model);
    reference.sites[0].quantum_numbers.erase("N");
    // Validate the complete reference even when the checked subset does not
    // contain its malformed site, or contains no operator factors at all.
    EXPECT_THROW(mdl::require_compatible_sites(reference, {{1, model.site_for(1)}}), std::invalid_argument);
    EXPECT_THROW(mdl::require_compatible_sites(reference, {}), std::invalid_argument);
    EXPECT_FALSE(reference.sites[0].quantum_numbers.count("N") != 0);
}
