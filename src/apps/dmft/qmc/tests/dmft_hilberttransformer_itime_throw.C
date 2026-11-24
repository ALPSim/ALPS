#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Behaviour lock for SemicircleHilbertTransformer::operator() on the
// imaginary-time (non-OMEGA_LOOP) self-consistency path.
//
// That imaginary-time Hilbert transform was never implemented: the body
// printed a message and called exit(1), followed by statically-unreachable
// code that default-constructed a std::shared_ptr<FourierTransformer> and
// immediately dereferenced it (a guaranteed null-deref had control ever
// reached it). The exit(1) is genuinely reachable — the Hirsch-Fye itime
// branch in main.C drives selfconsistency_loop, which calls this operator()
// — so a real DMFT run would hard-abort.
//
// This lock pins the post-fix contract: the unimplemented itime transform
// reports via std::logic_error (which the driver's catch surfaces cleanly)
// rather than terminating the process. It fails on the pre-fix exit(1)
// (the process exits non-zero before reaching the assertion) and passes
// once the body throws.

#include "hilberttransformer.h"

#include <alps/parameter.h>

#include <cstdio>
#include <stdexcept>

TEST(DmftRegression, hilberttransformer_itime_throw) {
  alps::Parameters parms;
  parms["FLAVORS"] = 2;
  parms["t"]       = 1.0;

  SemicircleHilbertTransformer transform(parms);

  // The values are irrelevant: the itime transform is unimplemented and
  // must report before touching them.
  itime_green_function_t G_tau(/*ntime*/ 11, /*nsite*/ 1, /*nflavor*/ 2);

  EXPECT_THROW(transform(G_tau, 0.0, 0.0, 10.0), std::logic_error);
}
