// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/include/alps/diag/product_problem.hpp.
// Adaptation: C++17 (no std::span), <alps/ed/...> include names, ALPS_ED_DECL
// exports, namespace alps::ed; observable assembly and expectation values added.
#pragma once

#include <alps/ed/eigensolver.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include <alps/ed_export.h>
#include <alps/models/hamiltonian.hpp>
#include <alps/models/physical_basis.hpp>

#include <complex>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace alps::ed {

struct ProductOptions {
    std::map<std::string, double> sector;
    std::size_t max_full_dimension = 10'000'000;
};
using EigenstateBuffer = std::variant<std::vector<double>, std::vector<std::complex<double>>>;
struct ProductEigenResult;

/// Immutable prepared Hermitian product problem. Copies share owned storage;
/// caller models/arrays can be changed or destroyed after preparation. Views
/// remain valid while a copy of the problem exists. Moved-from access throws.
/// Each solve owns its workspace and results and may run concurrently.
class ALPS_ED_DECL ProductProblem {
public:
    std::size_t dimension() const;
    std::vector<int> const& site_dims() const;
    /// Owned local labels, charges, parity and fermion convention.
    models::PhysicalBasis const& physical_basis() const;
    std::map<std::string, double> const& sector() const;
    /// Row-major product states, site 0 most significant, lexicographic order.
    std::vector<std::int32_t> const& states() const;
    /// Canonical CSR: matrix row is output, column is input. Exact zeros removed.
    std::vector<std::int64_t> const& row_offsets() const;
    std::vector<std::int64_t> const& column_indices() const;
    EigenstateBuffer const& matrix_values() const;
    /// Auto prefers Lanczos above dimension 1000. A requested iterative solve
    /// falls back to dense when the dimension does not exceed 8*count + 200,
    /// and an automatic SM request is always dense. No dense sparse-path temporary.
    ProductEigenResult eigenstates(EigenstateOptions const& = {}) const;
    /// Any operator over this problem's local spaces, assembled on its basis
    /// with the same parity strings. Transitions that leave the sector are
    /// dropped, so a charge-changing operator comes back projected.
    SparseMatrix<std::complex<double>> operator_matrix(models::CompiledOperator const&) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit ProductProblem(std::shared_ptr<const Data>);
    Data const& data() const;
    friend ALPS_ED_DECL ProductProblem prepare_product(models::TermSum const&,
        std::vector<models::SiteBasisPtr> const&, std::vector<std::string> const&,
        models::FermionConvention, ProductOptions const&);
};

struct ProductEigenResult {
    ProductProblem problem;
    std::vector<double> eigenvalues;
    /// Column-major; columns correspond to ascending eigenvalues.
    EigenstateBuffer eigenvectors;
    EigenstateMethod method = EigenstateMethod::dense;
    bool sparse_requested = false;
    std::optional<double> tolerance;
    std::optional<int> max_iterations, iterations;
    std::vector<double> residual_norms; // Lanczos only, in returned column order.
};

/// Validate the full unprojected operator's Hermiticity, then compile, enumerate
/// and assemble once. Requested charges must be declared; selecting any charge
/// validates all declarations globally. With no selected charge, unused false
/// conservation declarations do not prevent a full-space eigensolve. Locally
/// absent declared charges are neutral. No caller input is retained by reference.
ALPS_ED_DECL ProductProblem prepare_product(models::TermSum const&,
    std::vector<models::SiteBasisPtr> const&,
    std::vector<std::string> const& conserved = {},
    models::FermionConvention = models::FermionConvention::strict,
    ProductOptions const& = {});
ALPS_ED_DECL ProductProblem prepare_product(models::Hamiltonian const&, ProductOptions const& = {});

/// <psi_i|O|psi_i> for every eigenvector column of the result, with O
/// assembled by ProductProblem::operator_matrix.
ALPS_ED_DECL std::vector<std::complex<double>> expectation_values(
    ProductEigenResult const&, models::CompiledOperator const&);

} // namespace alps::ed
