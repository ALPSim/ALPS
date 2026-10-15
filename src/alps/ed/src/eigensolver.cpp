// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/src/fulldiag/eigensolver.cpp.
// Adaptation: C++17 (no std::span), <alps/ed/...> include names; LAPACK's
// divide-and-conquer drivers for every size (alpsx used Eigen below 128 and
// zheevr for complex matrices).
#include <alps/ed/eigensolver.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace alps::ed {
namespace {
// ALPS supports the LP64 LAPACK ABI. The trailing size_t arguments are the
// hidden character lengths modern Fortran compilers expect; implementations
// without them ignore the extra arguments.
extern "C" {
void dsyevd_(char const*, char const*, int const*, double*, int const*, double*,
             double*, int const*, int*, int const*, int*, std::size_t, std::size_t);
void zheevd_(char const*, char const*, int const*, std::complex<double>*, int const*, double*,
             std::complex<double>*, int const*, double*, int const*, int*, int const*, int*,
             std::size_t, std::size_t);
}

void check_status(int info, char const* routine) {
  if (info != 0)
    throw std::runtime_error(std::string(routine) + " failed with LAPACK status " + std::to_string(info));
}

int workspace_size(double size) {
  if (!std::isfinite(size) || size < 1. || size > std::numeric_limits<int>::max())
    throw std::length_error("dense eigensolver workspace exceeds the supported LAPACK integer range");
  return static_cast<int>(std::ceil(size));
}

template<class Scalar>
void lapack_solve(std::vector<Scalar>& matrix, int n, std::vector<double>& values, bool vectors) {
  char job = vectors ? 'V' : 'N', triangle = 'L';
  int info = 0, lwork = -1, liwork = -1, integer_query = 0;
  Scalar work_query{};
  if constexpr (std::is_same_v<Scalar, double>) {
    dsyevd_(&job, &triangle, &n, matrix.data(), &n, values.data(),
            &work_query, &lwork, &integer_query, &liwork, &info, 1, 1);
    check_status(info, "dsyevd workspace query");
    lwork = workspace_size(work_query);
    liwork = workspace_size(integer_query);
    std::vector<double> work(static_cast<std::size_t>(lwork));
    std::vector<int> integer_work(static_cast<std::size_t>(liwork));
    dsyevd_(&job, &triangle, &n, matrix.data(), &n, values.data(),
            work.data(), &lwork, integer_work.data(), &liwork, &info, 1, 1);
    check_status(info, "dsyevd");
  } else {
    int lrwork = -1;
    double real_query = 0.;
    zheevd_(&job, &triangle, &n, matrix.data(), &n, values.data(),
            &work_query, &lwork, &real_query, &lrwork, &integer_query, &liwork, &info, 1, 1);
    check_status(info, "zheevd workspace query");
    lwork = workspace_size(work_query.real());
    lrwork = workspace_size(real_query);
    liwork = workspace_size(integer_query);
    std::vector<Scalar> work(static_cast<std::size_t>(lwork));
    std::vector<double> real_work(static_cast<std::size_t>(lrwork));
    std::vector<int> integer_work(static_cast<std::size_t>(liwork));
    zheevd_(&job, &triangle, &n, matrix.data(), &n, values.data(),
            work.data(), &lwork, real_work.data(), &lrwork,
            integer_work.data(), &liwork, &info, 1, 1);
    check_status(info, "zheevd");
  }
}

template<class Scalar>
DenseEigenResult<Scalar> solve(std::vector<Scalar> matrix, std::size_t dimension, bool vectors,
                               std::optional<std::size_t> count, std::string const& which) {
  validate_eigensolver_controls(which);
  if (count && *count == 0) throw std::invalid_argument("eigenvalue count must be positive");
  if (dimension > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::length_error("dense eigensolver dimension exceeds the supported LAPACK integer range");
  if ((dimension == 0 && !matrix.empty()) ||
      (dimension != 0 && (matrix.size() / dimension != dimension || matrix.size() % dimension != 0)))
    throw std::invalid_argument("dense eigensolver matrix must contain dimension squared entries");
  for (auto value : matrix)
    if (!std::isfinite(std::real(value)) || !std::isfinite(std::imag(value)))
      throw std::invalid_argument("dense eigensolver matrix entries must be finite");

  DenseEigenResult<Scalar> result;
  result.eigenvalues.resize(dimension);
  if (dimension == 0) return result;
  // The divide-and-conquer workspace is counted in LP64 integers.
  auto const n64 = static_cast<std::int64_t>(dimension);
  if (dimension > static_cast<std::size_t>(std::numeric_limits<int>::max() / 26) ||
      (vectors && n64 * n64 > (std::numeric_limits<int>::max() - 1 - 6 * n64) / 2))
    throw std::length_error("dense eigensolver workspace exceeds the supported LAPACK integer range");
  lapack_solve(matrix, static_cast<int>(dimension), result.eigenvalues, vectors);
  if (vectors) result.eigenvectors_column_major = std::move(matrix);
  for (auto value : result.eigenvalues)
    if (!std::isfinite(value)) throw std::runtime_error("dense eigensolver returned non-finite eigenvalues");
  if (vectors) {
    for (auto value : result.eigenvectors_column_major)
      if (!std::isfinite(std::real(value)) || !std::isfinite(std::imag(value)))
        throw std::runtime_error("dense eigensolver returned non-finite eigenvectors");
  }
  if (count && *count < dimension) {
    const auto selected = select_eigenvalues(result.eigenvalues, *count, which);
    DenseEigenResult<Scalar> subset;
    for (auto index : selected) {
      subset.eigenvalues.push_back(result.eigenvalues[index]);
      if (vectors) {
        const auto first = result.eigenvectors_column_major.begin()
            + static_cast<std::ptrdiff_t>(index * dimension);
        subset.eigenvectors_column_major.insert(subset.eigenvectors_column_major.end(), first,
                                                first + static_cast<std::ptrdiff_t>(dimension));
      }
    }
    return subset;
  }
  return result;
}
} // namespace

void validate_eigensolver_controls(std::string const& which,
    std::optional<double> tolerance, std::optional<int> max_iterations) {
  if (which != "SA" && which != "LA" && which != "SM" && which != "LM")
    throw std::invalid_argument("which must be one of 'SA', 'LA', 'SM', 'LM'");
  if (tolerance && (!std::isfinite(*tolerance) || *tolerance <= 0.))
    throw std::invalid_argument("tolerance must be finite and positive");
  if (max_iterations && *max_iterations <= 0)
    throw std::invalid_argument("max_iterations must be positive");
}

std::vector<std::size_t> select_eigenvalues(std::vector<double> const& values,
    std::size_t count, std::string const& which) {
  validate_eigensolver_controls(which);
  if (count == 0) throw std::invalid_argument("eigenvalue count must be positive");
  for (double value : values)
    if (!std::isfinite(value)) throw std::invalid_argument("eigenvalues must be finite");
  std::vector<std::size_t> selected(values.size());
  std::iota(selected.begin(), selected.end(), 0);
  auto priority = [&](std::size_t a, std::size_t b) {
    if (which == "SA") return values[a] < values[b];
    if (which == "LA") return values[a] > values[b];
    if (std::abs(values[a]) == std::abs(values[b])) return values[a] < values[b];
    return which == "SM" ? std::abs(values[a]) < std::abs(values[b])
                         : std::abs(values[a]) > std::abs(values[b]);
  };
  // Stable sorts: ties keep the input order, so degenerate states come back
  // in the order the solver produced them.
  const auto size = std::min(count, values.size());
  std::stable_sort(selected.begin(), selected.end(), priority);
  selected.resize(size);
  std::stable_sort(selected.begin(), selected.end(), [&](auto a, auto b) { return values[a] < values[b]; });
  return selected;
}

DenseEigenResult<double> dense_eigensystem(
    std::vector<double> matrix, std::size_t dimension, bool eigenvectors,
    std::optional<std::size_t> count, std::string const& which) {
  return solve(std::move(matrix), dimension, eigenvectors, count, which);
}
DenseEigenResult<std::complex<double>> dense_eigensystem(
    std::vector<std::complex<double>> matrix, std::size_t dimension, bool eigenvectors,
    std::optional<std::size_t> count, std::string const& which) {
  return solve(std::move(matrix), dimension, eigenvectors, count, which);
}
} // namespace alps::ed
