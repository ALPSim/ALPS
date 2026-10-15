// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/tests/dense_eigensystem.cpp.
// Adaptation: <alps/ed/...> include names; converted to GoogleTest; residual
// and orthogonality checks written out instead of Eigen; the Eigen-vs-LAPACK
// benchmark is dropped.
#include <alps/ed/eigensolver.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace {
double const pi = std::acos(-1.);
double conjugate(double z) { return z; }
std::complex<double> conjugate(std::complex<double> z) { return std::conj(z); }

template<class T> void analytic() {
  for (int n : {0, 1, 2, 5, 32, 160}) {
    SCOPED_TRACE(n);
    auto const dim = static_cast<std::size_t>(n);
    std::vector<T> matrix(dim * dim, T{});   // column-major
    T hopping = -1.;
    if constexpr (!std::is_same_v<T, double>) hopping = {-1., 0.3};
    for (std::size_t i = 0; i < dim; ++i) {
      matrix[i * dim + i] = 0.125;
      if (i > 0) {
        matrix[(i - 1) * dim + i] = hopping;              // (i, i-1)
        if constexpr (std::is_same_v<T, double>) matrix[i * dim + i - 1] = hopping;
        else matrix[i * dim + i - 1] = std::conj(hopping); // (i-1, i)
      }
    }
    auto const original = matrix;
    for (bool vectors : {false, true}) {
      auto result = alps::ed::dense_eigensystem(matrix, dim, vectors);
      ASSERT_EQ(result.eigenvalues.size(), dim);
      ASSERT_EQ(result.eigenvectors_column_major.size(), vectors ? dim * dim : 0);
      EXPECT_EQ(matrix, original) << "dense eigensolver modified the input matrix";
      for (std::size_t j = 0; j < dim; ++j) {
        double expected = 0.125 - 2. * std::abs(hopping) * std::cos(static_cast<double>(j + 1) * pi / (n + 1));
        EXPECT_NEAR(result.eigenvalues[j], expected, 2e-12);
      }
      if (vectors && n != 0) {
        auto const& v = result.eigenvectors_column_major;
        for (std::size_t col = 0; col < dim; ++col) {
          for (std::size_t row = 0; row < dim; ++row) {
            T value = -result.eigenvalues[col] * v[col * dim + row];
            for (std::size_t k = 0; k < dim; ++k) value += matrix[k * dim + row] * v[col * dim + k];
            EXPECT_LE(std::abs(value), 2e-12);
          }
          for (std::size_t other = 0; other < dim; ++other) {
            T overlap{};
            for (std::size_t row = 0; row < dim; ++row) overlap += conjugate(v[col * dim + row]) * v[other * dim + row];
            EXPECT_LE(std::abs(overlap - T(col == other ? 1. : 0.)), 2e-12);
          }
        }
      }
    }
  }
}
} // namespace

TEST(DenseEigensystem, OpenChainSpectrumAndEigenvectorsRealAndComplex) {
  analytic<double>();
  analytic<std::complex<double>>();
}

TEST(DenseEigensystem, RejectsInvalidInput) {
  std::vector<double> wrong_size{1., 2., 3.};
  EXPECT_THROW(alps::ed::dense_eigensystem(wrong_size, 2, true), std::invalid_argument);
  EXPECT_THROW(alps::ed::dense_eigensystem(wrong_size, 0, true), std::invalid_argument);
  std::vector<double> nan{std::numeric_limits<double>::quiet_NaN()};
  EXPECT_THROW(alps::ed::dense_eigensystem(nan, 1, true), std::invalid_argument);
  std::vector<std::complex<double>> inf{{0., std::numeric_limits<double>::infinity()}};
  EXPECT_THROW(alps::ed::dense_eigensystem(inf, 1, true), std::invalid_argument);
  std::vector<double> one{1.};
  EXPECT_THROW(alps::ed::dense_eigensystem(one, 1, true, 0), std::invalid_argument);
  EXPECT_THROW(alps::ed::dense_eigensystem(one, 1, true, 1, "XX"), std::invalid_argument);
}

TEST(DenseEigensystem, SelectsRequestedSubsetInAscendingOrder) {
  std::vector<double> const diagonal{-9., -3., -.2, .7, 2., 6.};
  std::vector<double> matrix(36, 0.);
  for (std::size_t i = 0; i < 6; ++i) matrix[i * 6 + i] = diagonal[i];
  auto expect = [&](char const* which, std::vector<double> const& expected) {
    SCOPED_TRACE(which);
    auto const result = alps::ed::dense_eigensystem(matrix, 6, true, 2, which);
    ASSERT_EQ(result.eigenvalues.size(), 2);
    EXPECT_NEAR(result.eigenvalues[0], expected[0], 1e-14);
    EXPECT_NEAR(result.eigenvalues[1], expected[1], 1e-14);
    ASSERT_EQ(result.eigenvectors_column_major.size(), 12);
    for (std::size_t col = 0; col < 2; ++col) {
      double weight = 0.;
      for (std::size_t row = 0; row < 6; ++row)
        if (diagonal[row] == expected[col]) weight += std::abs(result.eigenvectors_column_major[col * 6 + row]);
      EXPECT_NEAR(weight, 1., 1e-14);
    }
  };
  expect("SA", {-9., -3.});
  expect("LA", {2., 6.});
  expect("SM", {-.2, .7});
  expect("LM", {-9., 6.});
  EXPECT_EQ(alps::ed::select_eigenvalues(diagonal, 10, "SA").size(), 6);
}
