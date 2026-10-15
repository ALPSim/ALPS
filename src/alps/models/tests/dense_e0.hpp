// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/tests/models/dense_e0.hpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports;
// dense reference uses LAPACK zheev instead of Eigen.
//
// Shared test helper: dense assembly of a CompiledModel emitted by
// alps::models::compile_hamiltonian(), and its ground-state energy. Per term
// the full-space operator is the Kronecker product over all sites of (factor
// op [right-multiplied by F when the site is on the term's JW string] | bare F
// on factor-less JW sites | identity); F = diag(parity_operator) of the
// site's basis. Site 0 is the most significant factor of the Kronecker order.
//
// The eigenvalues come from LAPACK's zheev; no Eigen dependency.

#ifndef ALPS_TESTS_MODELS_DENSE_E0_HPP
#define ALPS_TESTS_MODELS_DENSE_E0_HPP

#include <alps/models/hamiltonian.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

extern "C" void zheev_(const char* jobz, const char* uplo, const int* n,
                       std::complex<double>* a, const int* lda, double* w,
                       std::complex<double>* work, const int* lwork,
                       double* rwork, int* info);

namespace alps_test {

// Row-major square complex matrix.
struct DenseMatrix {
    std::size_t dim = 0;
    std::vector<std::complex<double>> data;
    explicit DenseMatrix(std::size_t n = 0) : dim(n), data(n * n) {}
    std::complex<double>& operator()(std::size_t r, std::size_t c) { return data[r * dim + c]; }
    std::complex<double> operator()(std::size_t r, std::size_t c) const { return data[r * dim + c]; }
    static DenseMatrix identity(std::size_t n) {
        DenseMatrix m(n);
        for (std::size_t i = 0; i < n; ++i) m(i, i) = 1.;
        return m;
    }
};

inline DenseMatrix kron(DenseMatrix const& a, DenseMatrix const& b) {
    DenseMatrix out(a.dim * b.dim);
    for (std::size_t r = 0; r < a.dim; ++r)
        for (std::size_t c = 0; c < a.dim; ++c) {
            auto const scale = a(r, c);
            if (scale == 0.) continue;
            for (std::size_t i = 0; i < b.dim; ++i)
                for (std::size_t j = 0; j < b.dim; ++j)
                    out(r * b.dim + i, c * b.dim + j) = scale * b(i, j);
        }
    return out;
}

inline double frobenius_distance(DenseMatrix const& a, DenseMatrix const& b) {
    if (a.dim != b.dim) throw std::invalid_argument("frobenius_distance: dimension mismatch");
    double sum = 0.;
    for (std::size_t i = 0; i < a.data.size(); ++i) sum += std::norm(a.data[i] - b.data[i]);
    return std::sqrt(sum);
}

inline DenseMatrix dense_hamiltonian(alps::models::Hamiltonian const& H,
                                     alps::models::CompiledModel const& d) {
    std::size_t const n = d.op.site_dims.size();
    std::size_t dim = 1;
    for (int sd : d.op.site_dims) dim *= static_cast<std::size_t>(sd);
    DenseMatrix total(dim);

    for (auto const& term : d.op.terms) {
        DenseMatrix acc = DenseMatrix::identity(1);
        for (std::size_t s = 0; s < n; ++s) {
            std::size_t const sd = static_cast<std::size_t>(d.op.site_dims[s]);
            DenseMatrix local = DenseMatrix::identity(sd);
            bool const in_jw = std::find(term.jordan_wigner.begin(), term.jordan_wigner.end(),
                                         static_cast<int>(s)) != term.jordan_wigner.end();
            auto const fit = std::find_if(term.factors.begin(), term.factors.end(),
                [s](auto const& factor) { return factor.site == static_cast<int>(s); });
            auto const& parity = H.site_for(static_cast<int>(s))->parity_operator;
            if (fit != term.factors.end()) {
                for (std::size_t r = 0; r < sd; ++r)
                    for (std::size_t c = 0; c < sd; ++c)
                        local(r, c) = fit->matrix[r * sd + c] * (in_jw ? parity[c] : 1.);
            } else if (in_jw) {
                for (std::size_t r = 0; r < sd; ++r) local(r, r) = parity[r];
            }
            acc = kron(acc, local);
        }
        for (std::size_t i = 0; i < total.data.size(); ++i)
            total.data[i] += term.coefficient * acc.data[i];
    }
    return total;
}

// Ascending eigenvalues of a Hermitian matrix.
inline std::vector<double> hermitian_eigenvalues(DenseMatrix matrix) {
    int const n = static_cast<int>(matrix.dim);
    std::vector<double> values(static_cast<std::size_t>(n));
    if (n == 0) return values;
    // zheev reads the lower triangle in column-major order, which is the
    // upper triangle of the row-major storage; the matrix is Hermitian, so
    // conjugate-transposing the storage gives LAPACK its expected layout.
    for (int r = 0; r < n; ++r)
        for (int c = r + 1; c < n; ++c)
            std::swap(matrix.data[static_cast<std::size_t>(r) * n + c],
                      matrix.data[static_cast<std::size_t>(c) * n + r]);
    for (auto& z : matrix.data) z = std::conj(z);
    int lwork = -1, info = 0;
    std::complex<double> query;
    std::vector<double> rwork(static_cast<std::size_t>(std::max(1, 3 * n - 2)));
    zheev_("N", "L", &n, matrix.data.data(), &n, values.data(), &query, &lwork, rwork.data(), &info);
    lwork = static_cast<int>(query.real());
    std::vector<std::complex<double>> work(static_cast<std::size_t>(std::max(1, lwork)));
    zheev_("N", "L", &n, matrix.data.data(), &n, values.data(), work.data(), &lwork, rwork.data(), &info);
    if (info != 0) throw std::runtime_error("zheev failed with info " + std::to_string(info));
    return values;
}

inline double ground_state_energy(alps::models::Hamiltonian const& H,
                                  alps::models::CompiledModel const& d) {
    return hermitian_eigenvalues(dense_hamiltonian(H, d)).front();
}

} // namespace alps_test

#endif // ALPS_TESTS_MODELS_DENSE_E0_HPP
