// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/include/alps/diag/projected_problem.hpp.
// Adaptation: C++17 (no std::span), <alps/ed/...> include names, ALPS_ED_DECL
// exports, namespace alps::ed; the finite point-group entry points
// (prepare_finite, finite_spectrum) are not ported.
#pragma once

#include <alps/ed/eigensolver.hpp>
#include <alps/ed/projected_basis.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include <alps/ed/spectrum.hpp>
#include <alps/ed_export.h>
#include <alps/models/hamiltonian.hpp>
#include <alps/models/physical_basis.hpp>

#include <complex>
#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace alps::ed {

struct MomentumSelection {
    // Row of the Cartesian momentum grid; the last direction varies fastest.
    std::size_t index = 0;
    // Canonical components 2*pi*j/period in [0, 2*pi), one per direction.
    std::vector<double> momentum;
};

// Snap a wavevector in radians onto the grid of translation periods. Each
// component is wrapped into [0, 2*pi) and must lie within 1e-10 radians of an
// allowed momentum 2*pi*j/period. Empty inputs select the trivial momentum.
// Periods above 2^52 or a grid exceeding size_t are rejected.
ALPS_ED_DECL MomentumSelection select_momentum(std::vector<std::size_t> const& periods,
                                               std::vector<double> const& wavevector);

struct ProjectedSpaces {
    // Declared conserved charges in sorted name order. Every local space has
    // one charge row per name, in this order.
    std::vector<std::string> charges;
    std::vector<SpectrumLocalSpace> local_spaces;
};

// Local parity and charge rows of a complete physical basis. Parity is the
// site parity under the strict fermion convention and +1 under the lattice
// convention, whose operators carry no exchange signs.
ALPS_ED_DECL ProjectedSpaces projected_spaces(models::PhysicalBasis const& physical);

struct ProjectedOptions {
    // Fixed totals of declared charges; omitted charges remain mixed.
    std::map<std::string, double> sector;
    // One wavevector component per translation, in radians. Leave it empty
    // without translations to select a charge-only product basis.
    std::vector<double> momentum;
};

struct ProjectedEigenResult;

/// Immutable prepared Hermitian problem in an owned projected basis.
/// Preparation validates the operator on the complete selected parent sector
/// and assembles the reduced matrix once; solves reuse it. Copies share owned
/// storage, no caller input is retained by reference and moved-from access
/// throws. Returned references remain valid while a copy of the problem exists.
/// Each solve owns its workspace and results and may run concurrently.
class ALPS_ED_DECL ProjectedProblem {
public:
    ProjectedBasis const& basis() const;
    std::size_t dimension() const;
    /// Names of the basis charge rows after model preparation; empty when an
    /// operator was prepared on an explicitly constructed basis.
    std::vector<std::string> const& charges() const;
    /// The prepared reduced Hamiltonian as canonical CSR (sorted columns,
    /// duplicates summed, exact zeros removed). Output states are rows; input
    /// states are columns.
    SparseMatrix<std::complex<double>> const& hamiltonian() const;
    /// Shares the model eigenstate policy of ProductProblem (see
    /// EigenstateOptions). Vectors are complex for every momentum.
    ProjectedEigenResult eigenstates(EigenstateOptions const& = {}) const;

private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit ProjectedProblem(std::shared_ptr<const Data>);
    Data const& data() const;
    static ProjectedProblem prepare(models::CompiledOperator const&, ProjectedBasis const&,
                                    std::vector<std::string> charges);
    friend ALPS_ED_DECL ProjectedProblem prepare_projected(models::CompiledOperator const&,
                                                           ProjectedBasis const&);
    friend ALPS_ED_DECL ProjectedProblem prepare_projected(models::CompiledOperator const&,
        ProjectedSpaces const&, std::vector<SpectrumTranslation> const&, ProjectedOptions const&);
    friend ALPS_ED_DECL ProjectedProblem prepare_projected(models::TermSum const&,
        std::vector<models::SiteBasisPtr> const&, std::vector<std::string> const&,
        models::FermionConvention, std::vector<SpectrumTranslation> const&, ProjectedOptions const&);
};

struct ProjectedEigenResult {
    // Retains the problem, including its reduced matrix.
    ProjectedProblem problem;
    // Ascending eigenvalues, actual method, block dimension, momentum and the
    // effective Lanczos controls; dense solves leave the controls empty.
    SpectrumBlock spectrum;
    // Columns in the problem basis, matching ascending eigenvalues.
    std::vector<std::complex<double>> eigenvectors_column_major;
    std::vector<double> residual_norms; // Lanczos only, in returned column order.
    bool sparse_requested = false;
};

// Validate an operator on an explicitly constructed basis: non-Hermiticity,
// false declared symmetries and nonfinite input or actions throw.
ALPS_ED_DECL ProjectedProblem prepare_projected(models::CompiledOperator const&,
                                                ProjectedBasis const&);

// Model-level compiled input: charge names must be sorted and unique, with one
// matching row in every local space. Validate Hermiticity and every declared
// charge on the full operator before selecting even an empty sector. Explicit
// Jordan-Wigner strings use local_spaces' parity. Translation invariance is
// checked on the selected parent sector. No full product matrix is constructed.
ALPS_ED_DECL ProjectedProblem prepare_projected(models::CompiledOperator const&,
    ProjectedSpaces const&, std::vector<SpectrumTranslation> const&,
    ProjectedOptions const& = {});

// Model preparation: validate the sector and momentum requests, then the full
// unprojected model's Hermiticity and conservation of every declared charge,
// because unselected charges remain mixed in the projected basis. Then build
// the basis from sorted charge rows and validate it as above. Translations
// map source to destination sites and must preserve the local spaces.
ALPS_ED_DECL ProjectedProblem prepare_projected(models::TermSum const&,
    std::vector<models::SiteBasisPtr> const&, std::vector<std::string> const& conserved,
    models::FermionConvention, std::vector<SpectrumTranslation> const& translations,
    ProjectedOptions const& = {});

// Translations are the model lattice's periodic-direction generators when a
// momentum is requested; an empty momentum selects a charge-only basis.
ALPS_ED_DECL ProjectedProblem prepare_projected(models::Hamiltonian const&,
                                                ProjectedOptions const& = {});

struct ModelSpectrumOptions {
    // Fixed totals of declared charges; every other declared charge is resolved.
    std::map<std::string, double> sector;
    // One wavevector component per translation selects that momentum;
    // all_momenta resolves the complete grid. Without either, no translations.
    std::vector<double> momentum;
    bool all_momenta = false;
    // Maximum levels per joint block; omitted for every level of every block.
    std::optional<int> count;
    std::string which = "SA";
    std::optional<double> tolerance;
    std::optional<int> max_iterations;
};

struct ModelSpectrum {
    // Sorted declared charge names, the order of SpectrumBlock::quantum_numbers.
    std::vector<std::string> charges;
    // Canonical selected momentum when one was requested; otherwise empty.
    std::vector<double> momentum;
    // Blocks follow symmetry_spectrum_from_local ordering, dispatch and defaults.
    std::vector<SpectrumBlock> blocks;
    // Resolved fixed totals and periodic direction labels. Raw translations use
    // positional labels; lattice adapters retain the actual lattice axes.
    std::map<std::string, double> sector;
    std::vector<int> directions;
};

// Joint charge/momentum spectrum of a model. Controls and the sector and
// momentum requests are checked first, then full-model Hermiticity and every
// declared charge, then translation invariance on the selected parent sector,
// as for prepare_projected. No parent product table or matrix is built.
ALPS_ED_DECL ModelSpectrum model_spectrum(models::TermSum const&,
    std::vector<models::SiteBasisPtr> const&, std::vector<std::string> const& conserved,
    models::FermionConvention, std::vector<SpectrumTranslation> const& translations,
    ModelSpectrumOptions const& = {});

// Translations are the model lattice's periodic-direction generators when a
// momentum or all momenta are requested.
ALPS_ED_DECL ModelSpectrum model_spectrum(models::Hamiltonian const&,
                                          ModelSpectrumOptions const& = {});

// Same full-model validation and block conventions as compiled preparation.
ALPS_ED_DECL ModelSpectrum model_spectrum(models::CompiledOperator const&,
    ProjectedSpaces const&, std::vector<SpectrumTranslation> const&,
    ModelSpectrumOptions const& = {});

} // namespace alps::ed
