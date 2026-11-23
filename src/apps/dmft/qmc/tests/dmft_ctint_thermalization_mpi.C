// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Under MPI the interaction-expansion solver stops at a check where the
// ranks' completed fractions sum to 1 or MAX_TIME has passed, and a rank that
// is still thermalizing then has no measurements. The solver must then throw
// on every rank: collecting results with only some ranks measured used to
// hang, and with no rank measured, rank 0 failed while the others waited.
//
// The first check comes after one update on each rank and the next one a
// minute later, so these tests end at the first check. To get there with one
// rank still thermalizing, they give the other ranks a far larger
// THERMALIZATION than rank 0, which is the state a slower rank reaches.

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <alps/solvers.hpp>

#include <mpi.h>

#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
const int n_meas = 25;

int rank() {
  int r = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &r);
  return r;
}

int size() {
  int s = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &s);
  return s;
}

// The solver appends to matrix_size and writes staggered_sz in the working
// directory.
struct working_directory {
  explicit working_directory(std::filesystem::path const & dir)
    : previous(std::filesystem::current_path()) { std::filesystem::current_path(dir); }
  ~working_directory() { std::filesystem::current_path(previous); }
  std::filesystem::path previous;
};

// Parameters for G0_f(i w_n) = 1/(i w_n + 0.4) on both flavors, written to
// this rank's directory.
alps::params ctint_parameters(std::filesystem::path const & dir, int sweeps, int thermalization,
                              int max_time = 60) {
  const double beta = 5.;
  const int n_matsubara = 32;
  const std::string g0 = (dir / "g0.h5").string();
  {
    alps::hdf5::archive ar(g0, "w");
    std::vector<std::complex<double> > g(n_matsubara);
    for (int n = 0; n < n_matsubara; ++n)
      g[n] = 1. / std::complex<double>(0.4, (2 * n + 1) * std::acos(-1.) / beta);
    ar << alps::make_pvp("/G0_0", g);
    ar << alps::make_pvp("/G0_1", g);
  }
  alps::params p;
  p.supply("INFILE", g0);
  p.supply("BETA", std::to_string(beta));
  p.supply("U", "4");
  p.supply("MU", "1");
  p.supply("ALPHA", "0.3");
  p.supply("N_ORBITALS", "2");
  p.supply("N_MATSUBARA", std::to_string(n_matsubara));
  p.supply("N_TAU", "64");
  p.supply("HISTOGRAM_MEASUREMENT", "1");
  p.supply("SWEEPS", std::to_string(sweeps));
  p.supply("THERMALIZATION", std::to_string(thermalization));
  p.supply("N_MEAS", std::to_string(n_meas));
  p.supply("MAX_TIME", std::to_string(max_time));
  p.supply("SEED", "3");
  return p;
}
}

TEST(DmftRegression, ctint_mpi_merges_thermalized_ranks) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  // Every rank thermalizes in its first update, and the fractions
  // (N_MEAS - THERMALIZATION) / SWEEPS sum to 1 at the first check.
  const int thermalization = 10;
  alps::params p = ctint_parameters(temporary.path(), size() * (n_meas - thermalization), thermalization);
  const std::string output = (temporary.path() / "out.h5").string();
  ASSERT_NO_THROW(alps::solvers::ctint(p, output));
  if (rank() == 0) {
    alps::hdf5::archive ar(output, "r");
    std::uint64_t count = 0;
    ar >> alps::make_pvp("/simulation/results/Sign/count", count);
    EXPECT_EQ(count, static_cast<std::uint64_t>(size()));
  }
}

TEST(DmftRegression, ctint_mpi_rank_still_thermalizing) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  // Rank 0 completes SWEEPS in its first update; the others are still
  // thermalizing at that check.
  alps::params p = ctint_parameters(temporary.path(), n_meas, rank() == 0 ? 0 : 2000000000);
  const std::string output = (temporary.path() / "out.h5").string();
  EXPECT_THROW(alps::solvers::ctint(p, output), std::runtime_error);
}

TEST(DmftRegression, ctint_mpi_no_rank_thermalized) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  // MAX_TIME=-1 makes the first check stop the run, with every rank still
  // thermalizing.
  alps::params p = ctint_parameters(temporary.path(), n_meas, 2000000000, -1);
  const std::string output = (temporary.path() / "out.h5").string();
  EXPECT_THROW(alps::solvers::ctint(p, output), std::runtime_error);
}
