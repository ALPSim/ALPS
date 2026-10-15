// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Measured observables of the classic DMRG solver against exact
// diagonalization: per-site and per-pair terms, expectation values of whole
// operators, factors that split into several charge components, scalar
// terms and several targeted states.
#include "support.hpp"
#include <gtest/gtest.h>

#include <cmath>

using namespace classic_test;

namespace {
void expect_rows(cd::Measurement const& measured, std::vector<double> const& reference, double tolerance,
                 std::size_t row = 0) {
    ASSERT_LT(row, measured.values.size());
    ASSERT_EQ(measured.values[row].size(), reference.size());
    for (std::size_t k = 0; k < reference.size(); ++k)
        EXPECT_NEAR(measured.values[row][k], reference[k], tolerance) << "column " << k;
}
} // namespace

TEST(ClassicDmrgMeasurements, SiteAndPairSpinTermsMatchExactDiagonalization) {
    auto const model = mdl::heisenberg(lat::chain(6, 1.0, lat::Boundary::open), 1.0);
    Sector const sector{{"Sz", 1.}};
    cd::Measurements requests;
    requests.terms["Sz"] = cd::site_terms(model, "Sz");
    requests.terms["SzSz"] = cd::pair_terms(model, "Sz", "Sz");
    requests.terms["SpSm"] = cd::pair_terms(model, "Sp", "Sm");
    requests.expectation["H"] = mdl::compile_hamiltonian(model).op;
    requests.expectation["Sz_total"] = cd::site_terms(model, "Sz");
    auto const result = cd::solve(model, exact_options(32, sector), requests);
    auto const states = ed_states(model, sector);
    ASSERT_EQ(result.measurements.size(), 5u);
    auto const& sz = result.measurements.at("Sz");
    EXPECT_EQ(sz.kind, cd::MeasurementKind::terms);
    ASSERT_EQ(sz.sites.size(), 6u);
    EXPECT_EQ(sz.sites[3], std::vector<int>{3});
    expect_rows(sz, ed_term_values(states, requests.terms["Sz"]), 1e-8);
    auto const& szsz = result.measurements.at("SzSz");
    ASSERT_EQ(szsz.sites.size(), 36u);
    EXPECT_EQ(szsz.sites[7], (std::vector<int>{1}));     // (1, 1) merges on site
    EXPECT_EQ(szsz.sites[8], (std::vector<int>{1, 2}));
    EXPECT_EQ(szsz.sites[13], (std::vector<int>{1, 2})); // (2, 1) in ascending site order
    expect_rows(szsz, ed_term_values(states, requests.terms["SzSz"]), 1e-8);
    expect_rows(result.measurements.at("SpSm"), ed_term_values(states, requests.terms["SpSm"]), 1e-8);
    auto const& energy = result.measurements.at("H");
    EXPECT_EQ(energy.kind, cd::MeasurementKind::expectation);
    ASSERT_EQ(energy.sites.size(), 1u);
    EXPECT_TRUE(energy.sites[0].empty());
    expect_rows(energy, {result.energies[0]}, 1e-8);
    expect_rows(result.measurements.at("Sz_total"), {1.}, 1e-8);
}

TEST(ClassicDmrgMeasurements, FermionDensitiesAndHoppingCorrelationsMatchExactDiagonalization) {
    auto const model = mdl::hubbard(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 4.0);
    Sector const sector{{"N", 4.}, {"Sz", 0.}};
    cd::Measurements requests;
    requests.terms["n"] = cd::site_terms(model, "n");
    requests.terms["nupndn"] = cd::site_terms(model, "n_updn");
    requests.terms["hop_up"] = cd::pair_terms(model, "c_dag_up", "c_up");
    requests.terms["nn"] = cd::pair_terms(model, "n", "n");
    requests.expectation["H"] = mdl::compile_hamiltonian(model).op;
    auto const result = cd::solve(model, exact_options(128, sector), requests);
    auto const states = ed_states(model, sector);
    expect_rows(result.measurements.at("n"), ed_term_values(states, requests.terms["n"]), 1e-7);
    expect_rows(result.measurements.at("nupndn"), ed_term_values(states, requests.terms["nupndn"]), 1e-7);
    expect_rows(result.measurements.at("hop_up"), ed_term_values(states, requests.terms["hop_up"]), 1e-7);
    expect_rows(result.measurements.at("nn"), ed_term_values(states, requests.terms["nn"]), 1e-7);
    expect_rows(result.measurements.at("H"), {result.energies[0]}, 1e-7);
}

TEST(ClassicDmrgMeasurements, FactorsSplitIntoChargeComponentsAndSectorLeavingPartsVanish) {
    // Sx Sx = (S+ S- + S- S+ + S+ S+ + S- S-) / 4: the last two leave the Sz
    // sector and contribute nothing, the first two are measured.
    auto const model = mdl::heisenberg(lat::chain(6, 1.0, lat::Boundary::open), 1.0);
    Sector const sector{{"Sz", 0.}};
    cd::Measurements requests;
    requests.terms["SxSx"] = cd::pair_terms(model, "Sx", "Sx");
    requests.expectation["SxSx_sum"] = cd::pair_terms(model, "Sx", "Sx");
    auto const result = cd::solve(model, exact_options(32, sector), requests);
    auto const states = ed_states(model, sector);
    expect_rows(result.measurements.at("SxSx"), ed_term_values(states, requests.terms["SxSx"]), 1e-8);
    expect_rows(result.measurements.at("SxSx_sum"), {ed_value(states, requests.terms["SxSx"])}, 1e-8);
}

TEST(ClassicDmrgMeasurements, ScalarTermsShiftEnergiesAndExpectationValues) {
    auto const site = mdl::spin(0.5);
    auto model = mdl::heisenberg(lat::chain(6, 1.0, lat::Boundary::open), 1.0);
    model.add(mdl::TermSum::located(site, "I", 0) * std::complex<double>(0.75));
    auto const plain = mdl::heisenberg(lat::chain(6, 1.0, lat::Boundary::open), 1.0);
    cd::Measurements requests;
    requests.expectation["shifted_identity"] =
        model.operator_decomposition(mdl::TermSum::located(site, "I", 2) * std::complex<double>(3.)).op;
    auto const result = cd::solve(model, exact_options(32, {{"Sz", 0.}}), requests);
    EXPECT_NEAR(result.energies[0], ed_energies(plain, {{"Sz", 0.}})[0] + 0.75, 1e-9);
    expect_rows(result.measurements.at("shifted_identity"), {3.}, 1e-12);
}

TEST(ClassicDmrgMeasurements, EveryTargetedStateGetsOneRow) {
    auto const model = mdl::heisenberg(lat::chain(6, 1.0, lat::Boundary::open), 1.0);
    Sector const sector{{"Sz", 0.}};
    cd::Measurements requests;
    requests.expectation["H"] = mdl::compile_hamiltonian(model).op;
    requests.terms["SzSz"] = cd::pair_terms(model, "Sz", "Sz");
    auto const result = cd::solve(model, exact_options(64, sector, 2), requests);
    auto const states = ed_states(model, sector, 2);
    ASSERT_EQ(result.energies.size(), 2u);
    ASSERT_EQ(result.measurements.at("H").values.size(), 2u);
    ASSERT_EQ(result.measurements.at("SzSz").values.size(), 2u);
    for (std::size_t row = 0; row < 2; ++row) {
        SCOPED_TRACE(row);
        expect_rows(result.measurements.at("H"), {result.energies[row]}, 1e-7, row);
        expect_rows(result.measurements.at("SzSz"), ed_term_values(states, requests.terms["SzSz"], row), 1e-7, row);
    }
}
