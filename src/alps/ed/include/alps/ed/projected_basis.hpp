// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/include/alps/diag/projected_basis.hpp.
// Adaptation: C++17 (no std::span), <alps/ed/...> include names, ALPS_ED_DECL
// exports, namespace alps::ed; the finite point-group bases
// (make_finite_group_basis, finite_spectrum_from_local) are not ported.
#pragma once

#include <alps/ed/product_basis.hpp>
#include <alps/ed/sparse_from_local.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include <alps/ed/spectrum.hpp>
#include <alps/ed_export.h>
#include <alps/models/compiled.hpp>

#include <complex>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace alps::ed {
namespace detail {
struct ProjectedData;
struct ProjectedAccess;
} // namespace detail

// An owned orthonormal symmetry basis. A row labels a reduced vector, not a
// product state. Copies share immutable storage; input arrays are copied.
// Translation eigenvalues are exp(+i*k). Site 0 is the most significant digit
// in the lexicographically ordered representatives. The returned references
// stay valid while a copy of the basis exists.
class ALPS_ED_DECL ProjectedBasis {
public:
    std::size_t dimension() const;
    std::size_t num_sites() const;
    std::size_t parent_dimension() const;
    std::size_t momentum_index() const;
    std::vector<int> const& site_dims() const;
    std::vector<SpectrumLocalSpace> const& local_spaces() const;
    std::vector<SpectrumTranslation> const& translations() const;
    std::vector<std::optional<double>> const& charge_targets() const;
    std::vector<double> const& momentum() const;
    /// Row-major [dimension, num_sites] orbit representatives.
    std::vector<std::int32_t> const& representatives() const;
    std::vector<std::size_t> const& orbit_sizes() const;

private:
    explicit ProjectedBasis(std::shared_ptr<const detail::ProjectedData> data);
    std::shared_ptr<const detail::ProjectedData> data_;
    friend struct detail::ProjectedAccess;
};

// Empty targets leave all charges unrestricted; otherwise provide one optional
// target per local-space charge row. Selected totals satisfy |Q-target|<1e-9.
// Empty translations select the product basis (momentum_index must be zero).
// The selected basis is built by streaming charge-accepted product states;
// neither the parent row table nor a parent Hamiltonian is materialized.
ALPS_ED_DECL ProjectedBasis make_projected_basis(
    std::vector<SpectrumLocalSpace> const& local_spaces,
    std::vector<SpectrumTranslation> const& translations,
    std::vector<std::optional<double>> const& charge_targets = {},
    std::size_t momentum_index = 0);

// Rectangular P_target^dagger O P_source. No symmetry or Hermiticity is
// required of O: local, momentum-changing and charge-changing observables are
// supported. Both bases must describe identical ordered physical local spaces.
// COO rows are target states, columns source states, and duplicate
// contributions sum.
ALPS_ED_DECL CooFromLocal projected_operator(models::CompiledOperator const& input,
    ProjectedBasis const& source, ProjectedBasis const& target);

/// Owned canonical CSR for the same rectangular projection. Scalar conversion
/// precedes reduction; no parent product matrix or row table is constructed.
template<class Scalar = std::complex<double>>
SparseMatrix<Scalar> projected_sparse_operator(models::CompiledOperator const& input,
    ProjectedBasis const& source, ProjectedBasis const& target) {
    return canonical_csr<Scalar>(projected_operator(input, source, target),
                                 target.dimension(), source.dimension());
}

// O applied to a state of the source basis, returned in the target basis.
// A column of ProjectedEigenResult::eigenvectors_column_major is passed as a
// copy of that column.
ALPS_ED_DECL std::vector<std::complex<double>> apply_projected_operator(
    models::CompiledOperator const& input, ProjectedBasis const& source,
    ProjectedBasis const& target, std::vector<std::complex<double>> const& state);
ALPS_ED_DECL std::complex<double> expect_projected_operator(
    models::CompiledOperator const& input, ProjectedBasis const& basis,
    std::vector<std::complex<double>> const& state);

struct ProjectedExpansion {
    ProductBasis basis;
    std::vector<std::complex<double>> amplitudes;
};
// Explicit opt-in expansion into the selected parent product sector. Ordinary
// projected operator measurements do not call this potentially large operation.
ALPS_ED_DECL ProjectedExpansion expand_projected_state(ProjectedBasis const& basis,
    std::vector<std::complex<double>> const& state);

// Joint charge/momentum spectra without parent CSR or parent product rows.
// Empty momentum_indices selects the full grid. Charge ordering/clustering and
// dense/Lanczos dispatch retain symmetry_spectrum's contract and defaults.
ALPS_ED_DECL std::vector<SpectrumBlock> symmetry_spectrum_from_local(
    models::CompiledOperator const& input,
    std::vector<SpectrumLocalSpace> const& local_spaces,
    std::vector<SpectrumTranslation> const& translations,
    std::vector<std::optional<double>> const& charge_targets, int count,
    std::vector<std::size_t> const& momentum_indices = {},
    std::string const& which = "SA", std::optional<double> tolerance = {},
    std::optional<int> max_iterations = {});

} // namespace alps::ed
