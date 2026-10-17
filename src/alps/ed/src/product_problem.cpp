// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/src/runtime/product_problem.cpp.
// Adaptation: C++17 (no std::span or templated lambdas), <alps/ed/...>
// include names, namespace alps::ed; observable assembly and expectation
// values added.
#include <alps/ed/product_problem.hpp>
#include <alps/ed/lanczos.hpp>
#include <alps/ed/product_basis.hpp>
#include <alps/ed/sparse_from_local.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include "eigenstate_policy.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace alps::ed {
struct ProductProblem::Data {
    std::vector<int> dims;
    models::PhysicalBasis physical;
    std::map<std::string, double> sector;
    ProductBasis basis;
    std::vector<std::int64_t> offsets, columns;
    EigenstateBuffer values;
    std::vector<std::vector<double>> parity;
};
ProductProblem::ProductProblem(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
ProductProblem::Data const& ProductProblem::data() const {
    if (!data_) throw std::logic_error("ProductProblem: moved-from problem");
    return *data_;
}
std::size_t ProductProblem::dimension() const { return data().basis.dimension; }
std::vector<int> const& ProductProblem::site_dims() const { return data().dims; }
models::PhysicalBasis const& ProductProblem::physical_basis() const { return data().physical; }
std::map<std::string, double> const& ProductProblem::sector() const { return data().sector; }
std::vector<std::int32_t> const& ProductProblem::states() const { return data().basis.states; }
std::vector<std::int64_t> const& ProductProblem::row_offsets() const { return data().offsets; }
std::vector<std::int64_t> const& ProductProblem::column_indices() const { return data().columns; }
EigenstateBuffer const& ProductProblem::matrix_values() const { return data().values; }

ProductProblem prepare_product(models::TermSum const& terms,
    std::vector<models::SiteBasisPtr> const& sites, std::vector<std::string> const& conserved,
    models::FermionConvention convention, ProductOptions const& options) {
    for (auto const& [name, target] : options.sector) {
        if (std::find(conserved.begin(), conserved.end(), name) == conserved.end())
            throw std::invalid_argument("sector contains a charge not declared conserved: " + name);
        if (!std::isfinite(target)) throw std::invalid_argument("sector totals must be finite");
    }
    auto physical = models::physical_basis(sites, conserved, convention);
    models::validate_hermiticity(terms, sites, convention);
    auto const compiled = models::compile_hamiltonian(terms, sites,
        options.sector.empty() ? std::vector<std::string>{} : conserved, convention);
    std::vector<BasisConstraint> constraints;
    for (auto const& [name, target] : options.sector) {
        auto& constraint = constraints.emplace_back();
        constraint.target = target;
        for (auto const& site : physical.sites) constraint.values.push_back(site.quantum_numbers.at(name));
    }
    auto result = std::make_shared<ProductProblem::Data>();
    result->physical = std::move(physical);
    result->sector = options.sector;
    result->dims = compiled.site_dims;
    result->basis = enumerate_product_basis(result->dims, constraints, options.max_full_dimension);
    if (result->basis.dimension > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::length_error("product problem exceeds the eigensolver dimension limit");
    for (auto const& site : sites) result->parity.push_back(site->parity_operator);
    std::vector<std::int64_t> const states(result->basis.states.begin(), result->basis.states.end());
    // Row by row: the uncombined contribution list of a sector is never held.
    std::visit([&](auto&& matrix) {
        result->offsets = std::move(matrix.offsets);
        result->columns = std::move(matrix.indices);
        result->values = std::move(matrix.values);
    }, build_canonical_csr(compiled, states, result->basis.dimension, result->parity));
    return ProductProblem(std::move(result));
}
ProductProblem prepare_product(models::Hamiltonian const& model, ProductOptions const& options) {
    std::vector<models::SiteBasisPtr> sites;
    for (int site = 0; site < model.num_sites(); ++site) sites.push_back(model.site_for(site));
    return prepare_product(model.expand(), sites, model.conserved_qns(), model.fermion_convention(), options);
}

ProductEigenResult ProductProblem::eigenstates(EigenstateOptions const& options) const {
    auto const& source = data();
    detail::validate_request(options);
    auto const n = dimension();
    auto const plan = detail::eigenstate_plan(options, n);
    ProductEigenResult result{*this, {}, {}, EigenstateMethod::dense, false, {}, {}, {}, {}};
    result.sparse_requested = detail::lanczos_requested(options.method, n);
    std::visit([&](auto const& values) {
        using Scalar = typename std::decay_t<decltype(values)>::value_type;
        if (plan.iterative) {
            auto computed = lanczos_csr(source.offsets, source.columns, values, static_cast<std::int64_t>(n),
                static_cast<int>(plan.count), options.which, plan.tolerance, plan.max_iterations, true);
            result.eigenvalues = std::move(computed.eigenvalues);
            result.eigenvectors = std::move(computed.eigenvectors);
            result.residual_norms = std::move(computed.residual_norms);
            result.method = EigenstateMethod::lanczos;
            result.tolerance = computed.tolerance;
            result.max_iterations = computed.max_iterations;
            result.iterations = computed.iterations;
        } else {
            std::vector<Scalar> dense;
            if (n && n > dense.max_size() / n) throw std::length_error("dense product matrix is too large");
            dense.resize(n * n);
            for (std::size_t row = 0; row < n; ++row)
                for (auto entry = source.offsets[row]; entry < source.offsets[row + 1]; ++entry)
                    dense[static_cast<std::size_t>(source.columns[static_cast<std::size_t>(entry)]) * n + row]
                        = values[static_cast<std::size_t>(entry)];
            // The dense API accepts a positive requested count even for an empty matrix.
            auto computed = dense_eigensystem(std::move(dense), n, true,
                                              static_cast<std::size_t>(options.count), options.which);
            result.eigenvalues = std::move(computed.eigenvalues);
            result.eigenvectors = std::move(computed.eigenvectors_column_major);
        }
    }, source.values);
    return result;
}

SparseMatrix<std::complex<double>> ProductProblem::operator_matrix(models::CompiledOperator const& op) const {
    auto const& source = data();
    if (op.site_dims != source.dims)
        throw std::invalid_argument("ProductProblem::operator_matrix: the operator's local dimensions differ from the problem's");
    std::vector<std::int64_t> const states(source.basis.states.begin(), source.basis.states.end());
    return std::get<SparseMatrix<std::complex<double>>>(
        build_canonical_csr(op, states, source.basis.dimension, source.parity, true));
}

std::vector<std::complex<double>> expectation_values(ProductEigenResult const& result,
                                                     models::CompiledOperator const& op) {
    auto const matrix = result.problem.operator_matrix(op);
    auto const n = result.problem.dimension();
    std::vector<std::complex<double>> out;
    std::visit([&](auto const& vectors) {
        if (n == 0) return;
        if (vectors.size() % n != 0)
            throw std::invalid_argument("expectation_values: eigenvector buffer does not match the problem dimension");
        auto const columns = vectors.size() / n;
        for (std::size_t column = 0; column < columns; ++column) {
            auto const* v = vectors.data() + column * n;
            std::complex<double> value{};
            for (std::size_t row = 0; row < n; ++row) {
                std::complex<double> product{};
                for (auto p = matrix.offsets[row]; p < matrix.offsets[row + 1]; ++p)
                    product += matrix.values[static_cast<std::size_t>(p)]
                        * std::complex<double>(v[static_cast<std::size_t>(matrix.indices[static_cast<std::size_t>(p)])]);
                value += std::conj(std::complex<double>(v[row])) * product;
            }
            out.push_back(value);
        }
    }, result.eigenvectors);
    return out;
}
} // namespace alps::ed
