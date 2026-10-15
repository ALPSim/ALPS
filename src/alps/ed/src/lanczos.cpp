// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Deflated Lanczos driver over the in-tree IETL (src/alps/ietl, the engine
// of the sparsediag application). Interface after skilledwolf/alpsx
// src/solvers/diag/src/sparsediag/lanczos_csr.ipp, whose Spectra
// implementation is not used here.
//
// One IETL run converges the lowest distinct eigenvalues of the operator and
// their Ritz vectors. Degenerate partners are invisible to a single Krylov
// sequence, so converged vectors are locked and the next run works on the
// deflated operator (csr_operator.hpp). Once `count` pairs are locked, one
// more run on the remaining space verifies that nothing lower was missed;
// a lower value replaces the largest locked one and the check repeats.

#include <alps/ed/lanczos.hpp>
#include <alps/ed/eigensolver.hpp>

#include "csr_operator.hpp"

#include <ietl/interface/ublas.h>
#include <ietl/vectorspace.h>
#include <ietl/lanczos.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

namespace alps::ed {
namespace {
using detail::CsrView;
using detail::DeflatedOperator;

template<class Scalar> using Vector = boost::numeric::ublas::vector<Scalar>;

struct Controls {
    double tolerance;
    int max_iterations;
};

// Deterministic start vectors: a given problem always produces the same
// Krylov sequence, and every run of the deflation loop gets a fresh one.
struct StartVector {
    std::mt19937_64 engine;
    explicit StartVector(std::uint64_t seed) : engine(seed) {}
    double operator()() { return std::uniform_real_distribution<double>(-1., 1.)(engine); }
};

template<class Scalar>
struct Harvest {
    std::vector<double> values;                      // converged, ascending, of the deflated operator
    std::vector<std::vector<Scalar>> vectors;        // matching Ritz vectors
    unsigned iterations = 0;
};

template<class Scalar>
bool finite_vector(std::vector<Scalar> const& v) {
    return std::all_of(v.begin(), v.end(), [](Scalar z) {
        return std::isfinite(std::real(z)) && std::isfinite(std::imag(z));
    });
}

// One IETL run: the `want` lowest distinct converged eigenvalues of `op` and
// their eigenvectors, as many as converged and were computed.
template<class Scalar>
Harvest<Scalar> run(DeflatedOperator<Scalar> const& op, std::int64_t dim, int want,
                    Controls const& controls, std::uint64_t seed) {
    using Space = ietl::vectorspace<Vector<Scalar>>;
    Space const space(static_cast<std::size_t>(dim));
    ietl::lanczos<DeflatedOperator<Scalar>, Space> lanczos(op, space);
    ietl::lanczos_iteration_nlowest<double> iteration(
        static_cast<unsigned>(controls.max_iterations), static_cast<unsigned>(want),
        controls.tolerance, controls.tolerance);
    StartVector const generator(seed);
    lanczos.calculate_eigenvalues(iteration, generator);

    Harvest<Scalar> out;
    out.iterations = iteration.iterations();
    if (lanczos.diagonal().size() < 2)
        throw std::runtime_error("lanczos_csr: did not converge: fewer than two Lanczos steps were taken");
    auto const& errors = lanczos.errors();
    auto const& values = lanczos.eigenvalues();
    for (std::size_t i = 0; i < values.size() && i < errors.size()
             && out.values.size() < static_cast<std::size_t>(want); ++i) {
        if (!std::isfinite(values[i]) || !std::isfinite(errors[i])) break;
        if (errors[i] > std::max(controls.tolerance, controls.tolerance * std::abs(values[i]))) break;
        out.values.push_back(values[i]);
    }
    if (out.values.empty())
        throw std::runtime_error("lanczos_csr: did not converge within "
            + std::to_string(controls.max_iterations) + " Lanczos steps"
            + (iteration.error_code() ? " (step limit reached)" : " (breakdown)"));

    // The targets are copied: IETL may extend its T matrix while computing
    // the vectors, which rebuilds the eigenvalue buffers.
    auto const targets = out.values;
    std::vector<Vector<Scalar>> vectors;
    ietl::Info<double> info;
    lanczos.eigenvectors(targets.begin(), targets.end(), std::back_inserter(vectors), info, generator);
    std::size_t computed = 0;
    for (; computed < targets.size() && computed < vectors.size(); ++computed) {
        if (info.error_info(static_cast<int>(computed)) != ietl::Info<double>::ok) break;
        std::vector<Scalar> vector(vectors[computed].begin(), vectors[computed].end());
        if (!finite_vector(vector)) break;
        out.vectors.push_back(std::move(vector));
    }
    out.values.resize(computed);
    if (computed == 0)
        throw std::runtime_error("lanczos_csr: IETL could not compute an eigenvector for the converged eigenvalue");
    return out;
}

// Orthonormalize v against the locked set (twice) and append it.
template<class Scalar>
void lock(std::vector<std::vector<Scalar>>& locked, std::vector<Scalar> v) {
    auto const n = v.size();
    for (int pass = 0; pass < 2; ++pass) detail::project_out(locked, v, n);
    double const norm = std::sqrt(std::real(detail::inner<Scalar>(v, v, n)));
    if (!(norm > 1e-3))
        throw std::runtime_error("lanczos_csr: a converged Ritz vector lies in the span of the locked eigenvectors");
    for (auto& z : v) z /= norm;
    locked.push_back(std::move(v));
}

template<class Scalar>
struct EndResult {
    std::vector<double> values;
    std::vector<std::vector<Scalar>> vectors;
    int iterations = 0;
    int runs = 0;
};

// The `count` lowest eigenpairs of sign*A, with multiplicity, in the order
// they were locked. `shift` lies above the spectrum of sign*A.
template<class Scalar>
EndResult<Scalar> lowest_eigenpairs(CsrView<Scalar> const& matrix, double sign, int count,
                                    Controls const& controls, double shift, double degenerate_tolerance) {
    std::vector<std::vector<Scalar>> locked;
    DeflatedOperator<Scalar> const op{&matrix, sign, shift, &locked, {}};
    EndResult<Scalar> out;
    std::uint64_t seed = 0x5EED0000u + (sign > 0 ? 0u : 0x10000u);
    while (true) {
        int const need = count - static_cast<int>(locked.size());
        auto const harvest = run(op, matrix.dim, std::max(need, 1), controls, seed++);
        ++out.runs;
        out.iterations += static_cast<int>(harvest.iterations);
        if (need > 0) {
            for (std::size_t i = 0; i < harvest.values.size(); ++i) {
                lock(locked, harvest.vectors[i]);
                out.values.push_back(harvest.values[i] + shift);
            }
            continue;
        }
        // Verification on the deflated space: a lower eigenvalue than the
        // largest locked one is a degenerate partner that replaces it.
        double const candidate = harvest.values.front() + shift;
        auto const worst = std::max_element(out.values.begin(), out.values.end());
        if (candidate < *worst - degenerate_tolerance) {
            auto const index = static_cast<std::size_t>(worst - out.values.begin());
            locked.erase(locked.begin() + static_cast<std::ptrdiff_t>(index));
            out.values.erase(worst);
            lock(locked, harvest.vectors.front());
            out.values.push_back(candidate);
            continue;
        }
        break;
    }
    out.vectors = std::move(locked);
    return out;
}

template<class Scalar>
std::pair<double, double> gershgorin_bounds(CsrView<Scalar> const& matrix) {
    double lower = std::numeric_limits<double>::infinity(), upper = -lower;
    for (std::int64_t row = 0; row < matrix.dim; ++row) {
        double diagonal = 0., radius = 0.;
        for (auto p = matrix.offsets[row]; p < matrix.offsets[row + 1]; ++p) {
            if (matrix.indices[p] == row) diagonal += detail::real_part(matrix.values[p]);
            else radius += std::abs(matrix.values[p]);
        }
        lower = std::min(lower, diagonal - radius);
        upper = std::max(upper, diagonal + radius);
    }
    return {lower, upper};
}

template<class Scalar>
bool is_diagonal(CsrView<Scalar> const& matrix) {
    for (std::int64_t row = 0; row < matrix.dim; ++row)
        for (auto p = matrix.offsets[row]; p < matrix.offsets[row + 1]; ++p)
            if (matrix.indices[p] != row) return false;
    return true;
}

template<class Scalar>
LanczosResult<Scalar> solve(std::vector<std::int64_t> const& offsets, std::vector<std::int64_t> const& indices,
                            std::vector<Scalar> const& values, std::int64_t dim, int count,
                            std::string const& which, double tolerance, int max_iterations, bool want_vectors) {
    if (which != "SA" && which != "LA" && which != "LM")
        throw std::invalid_argument(which == "SM"
            ? "lanczos_csr: smallest-magnitude eigenvalues are not available from Lanczos; use dense_eigensystem"
            : "lanczos_csr: which must be one of 'SA', 'LA', 'LM'");
    if (dim <= 0 || dim > std::numeric_limits<int>::max())
        throw std::invalid_argument("lanczos_csr: dimension is outside the supported range");
    if (count <= 0 || static_cast<std::int64_t>(count) >= dim)
        throw std::invalid_argument("lanczos_csr: require 0 < count < dimension; got count="
            + std::to_string(count) + ", dimension=" + std::to_string(dim));
    if (offsets.size() != static_cast<std::size_t>(dim) + 1)
        throw std::invalid_argument("lanczos_csr: offsets length must be dimension+1");
    if (indices.size() != values.size())
        throw std::invalid_argument("lanczos_csr: indices and values length mismatch");
    if (offsets.front() != 0 || offsets.back() != static_cast<std::int64_t>(values.size())
        || !std::is_sorted(offsets.begin(), offsets.end()))
        throw std::invalid_argument("lanczos_csr: offsets must be nondecreasing from zero to the values length");
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (indices[i] < 0 || indices[i] >= dim)
            throw std::invalid_argument("lanczos_csr: indices must lie within the matrix dimensions");
        if (!std::isfinite(std::real(values[i])) || !std::isfinite(std::imag(values[i])))
            throw std::invalid_argument("lanczos_csr: values must be finite");
    }
    if (!std::isfinite(tolerance))
        throw std::invalid_argument("lanczos_csr: tolerance must be finite");

    // IETL needs at least two steps before it can judge convergence.
    Controls const controls{tolerance > 0. ? tolerance : 1e-10,
        std::max(2, max_iterations > 0 ? max_iterations
                                       : static_cast<int>(std::min<std::int64_t>(10 * dim, 1000)))};
    LanczosResult<Scalar> out;
    out.tolerance = controls.tolerance;
    out.max_iterations = controls.max_iterations;
    CsrView<Scalar> const matrix{dim, offsets.data(), indices.data(), values.data()};
    auto const n = static_cast<std::size_t>(dim);

    std::vector<double> candidates;
    std::vector<std::vector<Scalar>> vectors;
    if (is_diagonal(matrix)) {
        // Exact eigenpairs; no Krylov sequence exists for the zero operator.
        candidates.assign(n, 0.);
        for (std::int64_t row = 0; row < dim; ++row)
            for (auto p = offsets[row]; p < offsets[row + 1]; ++p)
                candidates[static_cast<std::size_t>(row)] += detail::real_part(values[p]);
        auto const selected = select_eigenvalues(candidates, static_cast<std::size_t>(count), which);
        out.nconv = static_cast<int>(selected.size());
        for (auto index : selected) {
            out.eigenvalues.push_back(candidates[index]);
            if (want_vectors) {
                out.eigenvectors.insert(out.eigenvectors.end(), n, Scalar{});
                out.eigenvectors[out.eigenvectors.size() - n + index] = Scalar{1.};
                out.residual_norms.push_back(0.);
            }
        }
        return out;
    }

    auto const bounds = gershgorin_bounds(matrix);
    double const lower = bounds.first, upper = bounds.second;
    double const margin = 1. + 0.01 * (upper - lower);
    double const degenerate_tolerance = 100. * controls.tolerance * std::max(1., upper - lower);
    auto collect = [&](double sign) {
        double const shift = (sign > 0 ? upper : -lower) + margin;
        auto end = lowest_eigenpairs(matrix, sign, count, controls, shift, degenerate_tolerance);
        out.iterations += end.iterations;
        out.runs += end.runs;
        for (std::size_t i = 0; i < end.values.size(); ++i) {
            candidates.push_back(sign * end.values[i]);
            vectors.push_back(std::move(end.vectors[i]));
        }
    };
    if (which != "LA") collect(1.);
    if (which != "SA") collect(-1.);

    // Select `count` candidates by priority; a candidate inside the span of
    // the already selected ones is the same eigenpair found from both ends.
    std::vector<std::size_t> order(candidates.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        if (which == "SA") return candidates[a] < candidates[b];
        if (which == "LA") return candidates[a] > candidates[b];
        if (std::abs(candidates[a]) == std::abs(candidates[b])) return candidates[a] < candidates[b];
        return std::abs(candidates[a]) > std::abs(candidates[b]);
    });
    std::vector<std::vector<Scalar>> selected;
    std::vector<double> selected_values;
    for (auto index : order) {
        if (selected.size() == static_cast<std::size_t>(count)) break;
        auto vector = vectors[index];
        for (int pass = 0; pass < 2; ++pass) detail::project_out(selected, vector, n);
        double const norm = std::sqrt(std::real(detail::inner<Scalar>(vector, vector, n)));
        if (!(norm > 1e-3)) continue;
        for (auto& z : vector) z /= norm;
        selected.push_back(std::move(vector));
        selected_values.push_back(candidates[index]);
    }

    std::vector<std::size_t> ascending(selected.size());
    std::iota(ascending.begin(), ascending.end(), 0);
    std::stable_sort(ascending.begin(), ascending.end(),
                     [&](std::size_t a, std::size_t b) { return selected_values[a] < selected_values[b]; });
    out.nconv = static_cast<int>(ascending.size());
    std::vector<Scalar> product(n);
    for (auto index : ascending) {
        out.eigenvalues.push_back(selected_values[index]);
        if (!want_vectors) continue;
        auto const& vector = selected[index];
        matrix.apply(vector, product);
        double residual = 0.;
        for (std::size_t i = 0; i < n; ++i) residual += std::norm(product[i] - selected_values[index] * vector[i]);
        out.residual_norms.push_back(std::sqrt(residual));
        out.eigenvectors.insert(out.eigenvectors.end(), vector.begin(), vector.end());
    }
    return out;
}
} // namespace

LanczosResult<double> lanczos_csr(
    std::vector<std::int64_t> const& offsets, std::vector<std::int64_t> const& indices,
    std::vector<double> const& values, std::int64_t dimension, int count,
    std::string const& which, double tolerance, int max_iterations, bool want_vectors) {
    return solve(offsets, indices, values, dimension, count, which, tolerance, max_iterations, want_vectors);
}

LanczosResult<std::complex<double>> lanczos_csr(
    std::vector<std::int64_t> const& offsets, std::vector<std::int64_t> const& indices,
    std::vector<std::complex<double>> const& values, std::int64_t dimension, int count,
    std::string const& which, double tolerance, int max_iterations, bool want_vectors) {
    return solve(offsets, indices, values, dimension, count, which, tolerance, max_iterations, want_vectors);
}

} // namespace alps::ed
