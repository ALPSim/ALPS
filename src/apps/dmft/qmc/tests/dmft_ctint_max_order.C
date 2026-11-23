// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The interaction-expansion solver does not add a vertex at MAX_ORDER
// vertices, so it samples the series cut at that order (#389). That can be
// wanted, to compare with low-order methods, but it was silent. The fraction
// of the updates at MAX_ORDER is now the result "MaxOrderReached", and a
// nonzero fraction is reported with a warning on stderr.

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/params.hpp>
#include <alps/solvers.hpp>

#include <filesystem>
#include <string>

#ifdef ALPS_HAVE_MPI
#include <mpi.h>
#endif

namespace {
bool is_master() {
#ifdef ALPS_HAVE_MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank == 0;
#else
  return true;
#endif
}

struct working_directory {
  explicit working_directory(std::filesystem::path const & dir)
    : previous(std::filesystem::current_path()) { std::filesystem::current_path(dir); }
  ~working_directory() { std::filesystem::current_path(previous); }
  std::filesystem::path previous;
};

// The Hubbard atom at U = 2, beta = 4 has a mean order of 4 per flavor.
alps::params parameters(int max_order) {
  alps::params p;
  p.supply("ATOMIC", "1");
  p.supply("BETA", "4");
  p.supply("U", "2");
  p.supply("ALPHA", "-0.01");
  p.supply("N_TAU", "50");
  p.supply("N_MATSUBARA", "40");
  p.supply("SWEEPS", "100000");
  p.supply("THERMALIZATION", "1000");
  p.supply("MAX_TIME", "60");
  p.supply("SEED", "5");
  if (max_order > 0)
    p.supply("MAX_ORDER", std::to_string(max_order));
  return p;
}

struct run_result {
  double max_order_reached;
  std::string stderr_text;
};

run_result run(std::filesystem::path const & dir, int max_order) {
  const std::string output = (dir / ("max_order_" + std::to_string(max_order) + ".h5")).string();
  testing::internal::CaptureStderr();
  alps::solvers::ctint(parameters(max_order), output);
  run_result r;
  r.stderr_text = testing::internal::GetCapturedStderr();
  r.max_order_reached = -1;
  if (is_master()) {
    alps::hdf5::archive ar(output, "r");
    ar >> alps::make_pvp("/simulation/results/MaxOrderReached/mean/value", r.max_order_reached);
  }
  return r;
}
}

TEST(DmftRegression, ctint_max_order_reached) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  const run_result r = run(temporary.path(), 2);
  if (!is_master()) return;
  // about 0.57 of the updates at order 2
  EXPECT_GT(r.max_order_reached, 0.3);
  EXPECT_LT(r.max_order_reached, 0.8);
  EXPECT_NE(r.stderr_text.find("WARNING: CT-INT: the expansion order was at MAX_ORDER = 2 in "),
            std::string::npos) << r.stderr_text;
  EXPECT_NE(r.stderr_text.find("the results are those of the expansion up to order 2, not of the full model."),
            std::string::npos) << r.stderr_text;
}

TEST(DmftRegression, ctint_max_order_default) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  const run_result r = run(temporary.path(), 0);
  if (!is_master()) return;
  EXPECT_EQ(r.max_order_reached, 0.);
  EXPECT_EQ(r.stderr_text.find("MAX_ORDER"), std::string::npos) << r.stderr_text;
}
