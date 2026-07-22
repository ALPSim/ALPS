// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Serialization/behaviour lock for ExternalSolver::solve() (the
// imaginary-time path): the OUTFILE value persisted into the input
// archive's /parameters group must be the OUTPUT filename, not the input
// filename. solve() used to set p["OUTFILE"] = infile (a copy-paste slip),
// while its Matsubara twin solve_omega() correctly uses outfile.
//
// In-tree this is latent — the bundled solver exes take the output path
// from argv[2] and ignore the archive's OUTFILE — but an out-of-tree
// external solver that reads OUTFILE from the input archive would be
// handed the wrong (input) path. This test pins the persisted value.
//
// Mechanism: point the solver at a non-existent executable. solve() writes
// and closes the input archive (with /parameters), then call() runs
// std::system on the bogus command, gets a non-zero result and throws —
// before it would delete the input archive on a successful run. We then
// reopen the input archive and assert /parameters/OUTFILE == the output
// filename.

#include "externalsolver.h"
#include "green_function.h"

#include <alps/parameter.h>
#include <alps/hdf5.hpp>

#include <cstdio>
#include <boost/filesystem.hpp>
#include <stdexcept>
#include <string>

int main() {
  const std::string tmpbase =
      (boost::filesystem::temp_directory_path() /
       "alps_dmft_extsolver_outfile_param_lock").string();
  const std::string infile  = tmpbase + ".in.h5";
  const std::string outfile = tmpbase + ".out.h5";
  boost::system::error_code ec;
  boost::filesystem::remove(infile, ec);
  boost::filesystem::remove(outfile, ec);

  alps::Parameters parms;
  parms["TMPNAME"] = tmpbase;   // makes infile/outfile deterministic
  parms["N"]       = 10;        // n_tau
  parms["SITES"]   = 1;
  parms["FLAVORS"] = 2;
  // No SC_WRITE_DELTA: solve() writes /G0 directly (no bandstructure needed).

  itime_green_function_t G0(/*ntime*/ 11, /*nsite*/ 1, /*nflavor*/ 2);
  for (unsigned i = 0; i <= 10; ++i) { G0(i, 0) = -0.5; G0(i, 1) = -0.5; }

  // A solver path that does not exist -> std::system returns non-zero ->
  // call() throws after writing (but before deleting) the input archive.
  ExternalSolver solver(
      boost::filesystem::path("alps_nonexistent_solver_outfile_param_lock"));
  bool solve_threw = false;
  try {
    solver.solve(G0, parms);
  } catch (const std::exception&) {
    solve_threw = true;  // expected (bogus solver exe)
  }
  if (!solve_threw) {
    std::printf("FAIL: solve() unexpectedly succeeded with a bogus solver exe\n");
    return 1;
  }
  if (!boost::filesystem::exists(infile)) {
    std::printf("FAIL: input archive %s not preserved for inspection\n",
                infile.c_str());
    return 1;
  }

  alps::Parameters written;
  {
    alps::hdf5::archive ar(infile, "r");
    ar >> alps::make_pvp("/parameters", written);
  }
  boost::filesystem::remove(infile, ec);

  const std::string got = written["OUTFILE"];
  if (got != outfile) {
    std::printf("FAIL: persisted /parameters/OUTFILE = '%s', expected '%s' "
                "(itime solve() must set OUTFILE to the output filename, "
                "matching solve_omega())\n",
                got.c_str(), outfile.c_str());
    return 1;
  }
  std::printf("OK: persisted /parameters/OUTFILE = '%s'\n", got.c_str());
  return 0;
}
