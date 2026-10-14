/*
 * Copyright (C) 1998-2018 ALPS Collaboration.
 * All rights reserved. Use is subject to license terms. See ../LICENSE.TXT
 * For use in publications, see ../ACKNOWLEDGE.TXT
 */

// Grid checks adapted from CQMP/Maxent test/gridTest.cpp at
// a3992534228b0de732e033ae23154969724814be: combine the odd/even cases,
// supply ALPS parameters directly, and check finite, ordered grid points.
#include "maxent.hpp"
#include "eigen_lapack.hpp"
#include "maxent_grid.hpp"
#include "maxent_kernel.hpp"

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

TEST(MaxEntDefaultModel, TabulatedInterpolationAndExtrapolation) {
  alps::testing::TemporaryDirectory temporary;
  const auto filename = temporary.path() / "model.dat";
  {
    std::ofstream model(filename.string());
    model << "\n  # comment\n\n-2 0.08 trailing annotation\n0 0.10\n\n2 0.12\n\n";
  }
  alps::params parameters;
  parameters["DEFAULT_MODEL"] = filename.string();
  parameters["OMEGA_MIN"] = -3.0;
  parameters["OMEGA_MAX"] = 3.0;
  TabFunction model(parameters, "DEFAULT_MODEL");
  for (double omega : {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0})
    EXPECT_NEAR(model(omega), 0.10 + 0.01 * omega, 1e-14);
}

TEST(MaxEntDefaultModel, RejectInvalidTables) {
  alps::testing::TemporaryDirectory temporary;
  const auto filename = temporary.path() / "model.dat";
  alps::params parameters;
  parameters["DEFAULT_MODEL"] = filename.string();
  parameters["OMEGA_MIN"] = -2.0;
  parameters["OMEGA_MAX"] = 2.0;
  for (const char* contents : {"", "# comment only\n", "0 1\n",
                               "-2 1\ninvalid row\n2 1\n", "-2 1\n0\n2 1\n",
                               "-2 1\n-2 1\n", "2 1\n-2 1\n"}) {
    SCOPED_TRACE(contents);
    {
      std::ofstream model(filename.string());
      model << contents;
    }
    EXPECT_THROW(TabFunction(parameters, "DEFAULT_MODEL"), std::invalid_argument);
  }
}

TEST(MaxEntSVD, RectangularThinFactors) {
  matrix_type tall(5, 3);
  tall << 3, 0, 0,
          4, 0, 0,
          0, 0, 1,
          0, 3, 0,
          0, 0, 0;
  for (bool transpose : {false, true}) {
    const matrix_type original = transpose ? matrix_type(tall.transpose()) : tall;
    SCOPED_TRACE(transpose ? "wide" : "tall");
    matrix_type input = original;
    // Previous allocations assumed NDAT <= NFREQ, leaving singular values and
    // factor columns unwritten when the input kernel had more rows than columns.
    const auto rows = input.rows(), columns = input.cols();
    const double poison = std::numeric_limits<double>::quiet_NaN();
    vector_type singular = vector_type::Constant(rows, poison);
    matrix_type left = matrix_type::Constant(rows, rows, poison);
    matrix_type right = matrix_type::Constant(rows, columns, poison);
    lapack_svd(input, singular, right, left);
    ASSERT_EQ(singular.size(), 3);
    ASSERT_EQ(left.rows(), rows);
    ASSERT_EQ(left.cols(), 3);
    ASSERT_EQ(right.rows(), 3);
    ASSERT_EQ(right.cols(), columns);
    ASSERT_TRUE(singular.allFinite());
    ASSERT_TRUE(left.allFinite());
    ASSERT_TRUE(right.allFinite());
    EXPECT_NEAR(singular[0], 5.0, 1e-12);
    EXPECT_NEAR(singular[1], 3.0, 1e-12);
    EXPECT_NEAR(singular[2], 1.0, 1e-12);
    EXPECT_LT((left * singular.asDiagonal() * right - original).norm(), 1e-12);
  }
}

TEST(MaxEntGrid, EndpointsAndOrdering) {
  for (int count : {5, 31, 128, 129, 2000, 2001}) {
    for (const char* name : {"Lorentzian", "half lorentzian", "Quadratic", "log", "linear"}) {
      SCOPED_TRACE(std::string(name) + ": " + std::to_string(count));
      alps::params parameters;
      parameters["NFREQ"] = count;
      parameters["FREQUENCY_GRID"] = std::string(name);
      parameters["CUT"] = 0.01;
      parameters["SPREAD"] = 4.0;
      parameters["LOG_MIN"] = 1e-4;
      grid points(parameters);
      ASSERT_EQ(points.t_array().size(), static_cast<std::size_t>(count + 1));
      EXPECT_NEAR(points(0), 0.0, 1e-10);
      EXPECT_NEAR(points(count), 1.0, 1e-10);
      for (int index = 1; index <= count; ++index) {
        ASSERT_TRUE(std::isfinite(points(index)));
        ASSERT_GT(points(index), points(index - 1));
      }
      if (count % 2 == 0 && std::string(name) != "half lorentzian")
        EXPECT_NEAR(points(count / 2), 0.5, 1e-10);
    }
  }
}

TEST(MaxEntGrid, LogarithmicBoundsAndSymmetry) {
  alps::params parameters;
  parameters["FREQUENCY_GRID"] = std::string("log");
  parameters["CUT"] = 0.01;
  parameters["LOG_MIN"] = 1e-4;
  for (int count : {3, 4, 5, 31, 129, 2000, 2001}) {
    SCOPED_TRACE(count);
    parameters["NFREQ"] = count;
    grid points(parameters);
    ASSERT_EQ(points.t_array().size(), static_cast<std::size_t>(count + 1));
    EXPECT_NEAR(points(0), 0.0, 1e-14);
    EXPECT_NEAR(points(count), 1.0, 1e-14);
    for (int index=0; index<=count; ++index) {
      EXPECT_NEAR(points(index) + points(count-index), 1.0, 1e-14);
      if (index > 0)
        EXPECT_GT(points(index), points(index-1));
    }
    if (count % 2 != 0) {
      EXPECT_NEAR(points(count/2), 0.5 - 1e-4, 1e-14);
      EXPECT_NEAR(points(count/2+1), 0.5 + 1e-4, 1e-14);
    }
  }
  for (int count : {0, 1, 2}) {
    parameters["NFREQ"] = count;
    EXPECT_THROW(grid invalid(parameters), std::invalid_argument);
  }
  parameters["NFREQ"] = 5;
  for (double minimum : {-1.0, 0.0, 0.5, 1.0, std::numeric_limits<double>::quiet_NaN()}) {
    parameters["LOG_MIN"] = minimum;
    EXPECT_THROW(grid invalid(parameters), std::invalid_argument);
  }
}

// Exact zeroth Legendre moment and parity identities independently check the
// Legendre quadrature kernel, including its regular limit at zero frequency.
TEST(MaxEntKernel, FermionicLegendreMoments) {
  constexpr double beta = 2.5;
  alps::params parameters;
  parameters["NDAT"] = 4;
  parameters["NFREQ"] = 3;
  parameters["BETA"] = beta;
  parameters["DATASPACE"] = std::string("legendre");
  parameters["KERNEL"] = std::string("fermionic");
  parameters["PARTICLE_HOLE_SYMMETRY"] = false;
  vector_type frequencies(3);
  frequencies << -1.0, 0.0, 1.0;
  vector_type input_grid = vector_type::Zero(4);
  kernel continuation(parameters, frequencies, input_grid);
  ASSERT_EQ(continuation.getDataspaceType(), legendre_dataspace);
  ASSERT_EQ(continuation.getKernelType(), legendre_fermionic_kernel);
  const auto& values = continuation();
  ASSERT_EQ(values.rows(), 4);
  ASSERT_EQ(values.cols(), 3);
  ASSERT_TRUE(values.allFinite());
  EXPECT_NEAR(values(0, 0), -std::tanh(beta / 2.0), 1e-10);
  EXPECT_NEAR(values(0, 1), -beta / 2.0, 1e-10);
  EXPECT_NEAR(values(0, 2), -std::tanh(beta / 2.0), 1e-10);
  for (int order = 1; order < 4; ++order) {
    EXPECT_NEAR(values(order, 1), 0.0, 1e-10);
    EXPECT_NEAR(values(order, 0), (order % 2 ? -1.0 : 1.0) * values(order, 2), 1e-10);
  }
}

TEST(MaxEntKernel, ComplexBosonicConvention) {
  alps::params parameters;
  parameters["NDAT"] = 4;
  parameters["NFREQ"] = 2;
  parameters["BETA"] = 2.0 * std::acos(-1.0);
  parameters["DATASPACE"] = std::string("frequency");
  parameters["KERNEL"] = std::string("bosonic");
  parameters["PARTICLE_HOLE_SYMMETRY"] = false;
  vector_type frequencies(2);
  frequencies << -1.0, 2.0;
  vector_type input_grid = vector_type::Zero(4);
  kernel continuation(parameters, frequencies, input_grid);
  ASSERT_EQ(continuation.getKernelType(), frequency_bosonic_kernel);
  // Omega / (i*nu - Omega) at nu=0 and nu=1, with real/imaginary rows.
  matrix_type expected(4, 2);
  expected << -1.0, -1.0,
               0.0,  0.0,
              -0.5, -0.8,
               0.5, -0.4;
  ASSERT_TRUE(continuation().allFinite());
  EXPECT_LT((continuation() - expected).norm(), 1e-12);
}

TEST(MaxEntParameters, RejectUnsupportedContinuations) {
  alps::params parameters;
  parameters["NDAT"] = 4;
  parameters["NFREQ"] = 3;
  parameters["BETA"] = 2.5;
  parameters["DATASPACE"] = std::string("legendre");
  parameters["KERNEL"] = std::string("bosonic");
  parameters["PARTICLE_HOLE_SYMMETRY"] = false;
  vector_type frequencies(3);
  frequencies << -1.0, 0.0, 1.0;
  vector_type input_grid = vector_type::Zero(4);
  EXPECT_THROW(kernel(parameters, frequencies, input_grid), std::invalid_argument);

  parameters["BACKCONTINUE"] = true;
  parameters["COVARIANCE_MATRIX"] = std::string("HDF5");
  EXPECT_THROW(MaxEntSimulation::define_parameters(parameters), std::invalid_argument);

  parameters["DATA_IN_HDF5"] = true;
  parameters["COVARIANCE_MATRIX"] = std::string("");
  EXPECT_THROW(MaxEntSimulation::define_parameters(parameters), std::invalid_argument);
}
