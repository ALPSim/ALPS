// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The catalog importer plus ALPS::ed against the legacy hamiltonian_matrix on
// the same catalog inputs and the same legacy basis enumeration: the matrices
// must agree element by element (fermionic signs included) and so must the
// complete spectra. The inputs follow the diag applications' sample
// parameter files (src/apps/test/parm*.in), sector by sector.
#include <alps/ed/eigensolver.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include <alps/model/catalog.hpp>
#include <alps/model/hamiltonian_matrix.hpp>
#include <alps/numeric/matrix.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace {
using Params = std::vector<std::pair<std::string, std::string>>;
struct Case {
    std::string name;
    Params params;
};
std::vector<Case> const& cases() {
    static std::vector<Case> const r = {
        {"spin chain, full space",
         {{"MODEL", "spin"}, {"LATTICE", "chain lattice"}, {"L", "4"}, {"J", "1"}, {"Jxy", "0.5"},
          {"CONSERVED_QUANTUMNUMBERS", "Sz"}}},
        {"spin chain, Sz=0",
         {{"MODEL", "spin"}, {"LATTICE", "chain lattice"}, {"L", "6"}, {"J", "1"}, {"Jxy", "0.5"},
          {"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"Sz_total", "0"}}},
        {"spin-1 ladder with anisotropy and field, Sz=1",
         {{"MODEL", "spin"}, {"LATTICE", "ladder"}, {"L", "3"}, {"local_S", "1"}, {"J0", "1"}, {"J1", "0.4"},
          {"D", "0.2"}, {"h", "0.1"}, {"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"Sz_total", "1"}}},
        {"hardcore bosons, N=2",
         {{"MODEL", "hardcore boson"}, {"LATTICE", "chain lattice"}, {"L", "4"}, {"t", "1"}, {"V", "0.5"},
          {"CONSERVED_QUANTUMNUMBERS", "N"}, {"N_total", "2"}}},
        {"boson Hubbard, N=3",
         {{"MODEL", "boson Hubbard"}, {"LATTICE", "chain lattice"}, {"L", "4"}, {"Nmax", "2"}, {"t", "1"},
          {"U", "4"}, {"CONSERVED_QUANTUMNUMBERS", "N"}, {"N_total", "3"}}},
        {"spinless fermions, periodic, N=2",
         {{"MODEL", "spinless fermions"}, {"LATTICE", "chain lattice"}, {"L", "4"}, {"t", "1"}, {"V", "0.5"},
          {"CONSERVED_QUANTUMNUMBERS", "N"}, {"N_total", "2"}}},
        {"spinless fermions, open, N=3",
         {{"MODEL", "spinless fermions"}, {"LATTICE", "open chain lattice"}, {"L", "5"}, {"t", "1"},
          {"V", "0.5"}, {"CONSERVED_QUANTUMNUMBERS", "N"}, {"N_total", "3"}}},
        {"fermion Hubbard chain, Nup=2 Ndown=1",
         {{"MODEL", "fermion Hubbard"}, {"LATTICE", "chain lattice"}, {"L", "4"}, {"t", "1"}, {"U", "4"},
          {"CONSERVED_QUANTUMNUMBERS", "Nup,Ndown"}, {"Nup_total", "2"}, {"Ndown_total", "1"}}},
        {"fermion Hubbard square, Nup=2 Ndown=2",
         {{"MODEL", "fermion Hubbard"}, {"LATTICE", "square lattice"}, {"L", "2"}, {"W", "2"}, {"t", "1"},
          {"U", "4"}, {"CONSERVED_QUANTUMNUMBERS", "Nup,Ndown"}, {"Nup_total", "2"}, {"Ndown_total", "2"}}},
        {"Ian's fermion Hubbard, N=2 Sz=0",
         {{"MODEL", "Ian's fermion Hubbard"}, {"LATTICE", "chain lattice"}, {"L", "6"}, {"t", "1"}, {"U", "2"},
          {"CONSERVED_QUANTUMNUMBERS", "N,Sz"}, {"N_total", "2"}, {"Sz_total", "0"}}},
        {"t-J chain, N=3 Sz=1/2",
         {{"MODEL", "t-J"}, {"LATTICE", "chain lattice"}, {"L", "4"}, {"t", "1"}, {"J", "0.4"},
          {"CONSERVED_QUANTUMNUMBERS", "N,Sz"}, {"N_total", "3"}, {"Sz_total", "0.5"}}},
    };
    return r;
}

std::vector<double> sorted_spectrum(std::vector<double> matrix, std::size_t dimension) {
    return alps::ed::dense_eigensystem(std::move(matrix), dimension, false).eigenvalues;
}
} // namespace

TEST(EdVsHamiltonianMatrix, CompiledCatalogOperatorsReproduceTheLegacyMatrices) {
    for (auto const& item : cases()) {
        SCOPED_TRACE(item.name);
        alps::Parameters p;
        for (auto const& [key, value] : item.params) p[key] = value;

        alps::hamiltonian_matrix<alps::numeric::matrix<double>> legacy(p);
        auto const& reference = legacy.matrix();
        auto const dimension = legacy.dimension();
        ASSERT_GT(dimension, 1u);
        auto const& basis = legacy.states_vector();
        ASSERT_EQ(basis.size(), dimension);
        std::size_t const sites = basis[0].size();
        std::vector<std::int64_t> states;
        for (auto const& state : basis)
            for (auto digit : state) states.push_back(digit);

        alps::models::CatalogModel<> model(p);
        ASSERT_EQ(model.basis().size(), sites);
        auto const op = model.compile_hamiltonian();
        auto const matrix = alps::ed::build_csr_from_local<std::complex<double>>(
            op, states, dimension, model.local_parity());
        std::vector<double> dense(dimension * dimension), expected(dimension * dimension);
        for (std::size_t row = 0; row < dimension; ++row)
            for (auto i = matrix.offsets[row]; i < matrix.offsets[row + 1]; ++i) {
                auto const value = matrix.values[static_cast<std::size_t>(i)];
                EXPECT_EQ(value.imag(), 0.);
                dense[static_cast<std::size_t>(matrix.indices[static_cast<std::size_t>(i)]) * dimension + row] =
                    value.real();
            }
        double scale = 1.;
        for (std::size_t row = 0; row < dimension; ++row)
            for (std::size_t column = 0; column < dimension; ++column) {
                expected[column * dimension + row] = reference(row, column);
                scale = std::max(scale, std::abs(reference(row, column)));
            }
        double difference = 0., asymmetry = 0.;
        for (std::size_t row = 0; row < dimension; ++row)
            for (std::size_t column = 0; column < dimension; ++column) {
                difference = std::max(difference,
                    std::abs(dense[column * dimension + row] - expected[column * dimension + row]));
                asymmetry = std::max(asymmetry,
                    std::abs(expected[column * dimension + row] - expected[row * dimension + column]));
            }
        EXPECT_LE(asymmetry, 1e-12 * scale);
        EXPECT_LE(difference, 1e-12 * scale);

        auto const spectrum = sorted_spectrum(dense, dimension), reference_spectrum = sorted_spectrum(expected, dimension);
        ASSERT_EQ(spectrum.size(), reference_spectrum.size());
        for (std::size_t i = 0; i < spectrum.size(); ++i)
            EXPECT_NEAR(spectrum[i], reference_spectrum[i], 1e-10 * scale);
    }
}
