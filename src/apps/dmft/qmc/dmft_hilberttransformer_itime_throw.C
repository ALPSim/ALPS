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

int main() {
  alps::Parameters parms;
  parms["FLAVORS"] = 2;
  parms["t"]       = 1.0;

  SemicircleHilbertTransformer transform(parms);

  // The values are irrelevant: the itime transform is unimplemented and
  // must report before touching them.
  itime_green_function_t G_tau(/*ntime*/ 11, /*nsite*/ 1, /*nflavor*/ 2);

  bool threw_logic_error = false;
  try {
    (void)transform(G_tau, /*mu*/ 0.0, /*h*/ 0.0, /*beta*/ 10.0);
  } catch (const std::logic_error&) {
    threw_logic_error = true;
  } catch (...) {
    std::printf("FAIL: operator() threw a non-logic_error exception\n");
    return 1;
  }

  if (!threw_logic_error) {
    std::printf("FAIL: SemicircleHilbertTransformer::operator() did not throw "
                "std::logic_error for the unimplemented itime path\n");
    return 1;
  }
  std::printf("OK: unimplemented itime Hilbert transform reports via "
              "std::logic_error\n");
  return 0;
}
