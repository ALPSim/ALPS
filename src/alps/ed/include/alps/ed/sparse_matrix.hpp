// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/include/alps/diag/sparse_matrix.hpp.
// Adaptation: C++17 (no std::span), <alps/ed/...> include names.
#pragma once

#include <alps/ed/sparse_from_local.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace alps::ed {

/// Owned canonical CSR. Rows are output states, columns input states. Columns
/// are sorted within each row, duplicates are summed, and exact zeros removed.
/// Scalar conversion precedes reduction, including for float/complex<float>.
/// Buffers are ordinary values: copies and later input edits are independent.
template<class Scalar> struct SparseMatrix {
    std::size_t rows = 0, columns = 0;
    std::vector<std::int64_t> offsets, indices;
    std::vector<Scalar> values;
};

namespace detail {
template<class Scalar> bool finite_scalar(Scalar value) {
    if constexpr (std::is_arithmetic_v<Scalar>) return std::isfinite(value);
    else return std::isfinite(value.real()) && std::isfinite(value.imag());
}
template<class Scalar> Scalar matrix_cast(std::complex<double> value) {
    if (!finite_scalar(value)) throw std::invalid_argument("sparse matrix contribution must be finite");
    if constexpr (std::is_arithmetic_v<Scalar>) {
        if (value.imag() != 0.)
            throw std::invalid_argument("complex matrix element cannot be cast to a real dtype: nonzero imaginary part");
        if constexpr (std::is_same_v<Scalar, bool>) return value.real() != 0.;
        else if constexpr (std::is_integral_v<Scalar>) {
            // Truncate then wrap modulo the scalar width. Avoid undefined
            // float-to-integer conversions and signed accumulation overflow.
            using Unsigned = std::make_unsigned_t<Scalar>;
            double magnitude = std::fmod(std::abs(std::trunc(value.real())),
                std::ldexp(1., std::numeric_limits<Unsigned>::digits));
            auto bits = static_cast<Unsigned>(magnitude);
            if (value.real() < 0.) bits = Unsigned(0) - bits;
            return static_cast<Scalar>(bits);
        } else return static_cast<Scalar>(value.real());
    } else return Scalar(value.real(), value.imag());
}
template<class Scalar> Scalar matrix_add(Scalar a, Scalar b) {
    if constexpr (std::is_same_v<Scalar, bool>) return a || b;
    else if constexpr (std::is_integral_v<Scalar>) {
        using Unsigned = std::make_unsigned_t<Scalar>;
        return static_cast<Scalar>(static_cast<Unsigned>(a) + static_cast<Unsigned>(b));
    } else return a + b;
}
} // namespace detail

/// Combine rectangular COO contributions in row-bucket/column-sort order.
/// Duplicates of one entry are summed in their COO order (the sort is stable).
/// Reject malformed indices and nonfinite input or floating-point results.
/// Real scalars require every contribution's imaginary part to be exactly zero,
/// even if those parts would cancel. Integer conversion truncates and wraps
/// modulo the scalar width; integer sums wrap and bool sums use logical OR.
/// The input's is_complex hint does not override the explicitly selected Scalar.
template<class Scalar>
SparseMatrix<Scalar> canonical_csr(CooFromLocal const& coo, std::size_t rows, std::size_t columns) {
    auto const limit = static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max());
    if (rows >= limit || columns > limit || coo.data.size() > limit)
        throw std::length_error("sparse matrix exceeds its index range");
    if (coo.rows.size() != coo.data.size() || coo.cols.size() != coo.data.size())
        throw std::invalid_argument("sparse matrix COO buffer lengths differ");
    SparseMatrix<Scalar> result{rows, columns, {}, {}, {}};
    result.offsets.assign(rows + 1, 0);
    for (std::size_t i = 0; i < coo.data.size(); ++i) {
        if (coo.rows[i] < 0 || static_cast<std::size_t>(coo.rows[i]) >= rows ||
            coo.cols[i] < 0 || static_cast<std::size_t>(coo.cols[i]) >= columns)
            throw std::invalid_argument("sparse matrix COO index outside shape");
        ++result.offsets[static_cast<std::size_t>(coo.rows[i]) + 1];
    }
    std::partial_sum(result.offsets.begin(), result.offsets.end(), result.offsets.begin());
    auto next = result.offsets;
    std::vector<std::pair<std::int64_t, Scalar>> entries(coo.data.size());
    for (std::size_t i = 0; i < entries.size(); ++i) {
        auto value = detail::matrix_cast<Scalar>(coo.data[i]);
        if (!detail::finite_scalar(value)) throw std::invalid_argument("sparse matrix scalar conversion overflowed");
        entries[static_cast<std::size_t>(next[static_cast<std::size_t>(coo.rows[i])]++)] = {coo.cols[i], value};
    }
    for (std::size_t row = 0; row < rows; ++row) {
        auto first = entries.begin() + result.offsets[row], last = entries.begin() + result.offsets[row + 1];
        std::stable_sort(first, last, [](auto const& a, auto const& b) { return a.first < b.first; });
        next[row] = static_cast<std::int64_t>(result.values.size());
        while (first != last) {
            auto column = first->first;
            Scalar value = first++->second;
            while (first != last && first->first == column) value = detail::matrix_add(value, first++->second);
            if (!detail::finite_scalar(value)) throw std::invalid_argument("sparse matrix accumulation overflowed");
            if (value != Scalar{}) { result.indices.push_back(column); result.values.push_back(value); }
        }
    }
    next[rows] = static_cast<std::int64_t>(result.values.size());
    result.offsets = std::move(next);
    return result;
}

/// Assemble an arbitrary operator on the caller's ordered, distinct product
/// rows. Non-Hermitian and charge-changing operators are supported. This owns
/// reduction; local-action and basis validation are shared with COO assembly.
template<class Scalar = std::complex<double>>
SparseMatrix<Scalar> build_csr_from_local(models::CompiledOperator const& input,
    std::vector<std::int64_t> const& states, std::size_t count,
    std::vector<std::vector<double>> const& parity = {}) {
    return canonical_csr<Scalar>(build_coo_from_local(input, states, count, parity), count, count);
}

using CanonicalCsr = std::variant<SparseMatrix<double>, SparseMatrix<std::complex<double>>>;

/// The canonical CSR of build_coo_from_local's contributions, assembled one
/// output row at a time: offsets, indices and values equal
/// canonical_csr<S>(build_coo_from_local(input, states, count, parity)), with S
/// real exactly when the COO's is_complex flag is false (or complex when
/// `complex_values`). Each row comes from the operator's adjoint applied to
/// that row's state (conjugated coefficients, conjugate-transposed factors,
/// factor-site string signs folded in; all exact), so the uncombined
/// contribution list is never held: peak memory is the result plus one row.
/// Same validation and exceptions as the COO path.
ALPS_ED_DECL CanonicalCsr build_canonical_csr(models::CompiledOperator const& input,
    std::vector<std::int64_t> const& states, std::size_t count,
    std::vector<std::vector<double>> const& parity = {}, bool complex_values = false);

} // namespace alps::ed
