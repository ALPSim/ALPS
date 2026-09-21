// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Lock for the dmft driver's CT-INT (Interaction Expansion) run-type
// selection on the Matsubara (OMEGA_LOOP) path.
//
// main.C used to select the factory with a tangled if / standalone-if /
// else-if chain that (a) built the single-site half-filling solver for
// FLAVORS==1 and then immediately overwrote it (a discarded construction),
// (b) compared FLAVORS as a string, and (c) fell through to
// ExternalSolver("Interaction Expansion") — a guaranteed exec failure — for
// any other (FLAVORS, SITES) combination. The decision is now a single
// mutually-exclusive function, select_interaction_expansion(flavors, sites),
// that this test pins. The mapping preserves the previous *effective*
// selection (FLAVORS==2 -> single site; else SITES==1 -> multiband) and turns
// the unreachable-as-a-real-solver combination into an explicit `unsupported`.

#include "interaction_expansion_choice.h"

#include <cstdio>

int main() {
  struct Case {
    int flavors;
    int sites;
    interaction_expansion_choice expect;
    const char* desc;
  };
  const Case cases[] = {
    {2, 1, interaction_expansion_choice::single_site_hubbard, "FLAVORS=2,SITES=1 -> single site"},
    {2, 4, interaction_expansion_choice::single_site_hubbard, "FLAVORS=2 ignores SITES"},
    {1, 1, interaction_expansion_choice::multiband_density,   "FLAVORS=1,SITES=1 -> multiband (was sshf-then-mbd)"},
    {4, 1, interaction_expansion_choice::multiband_density,   "FLAVORS=4,SITES=1 -> multiband"},
    {1, 2, interaction_expansion_choice::unsupported,         "FLAVORS=1,SITES=2 -> unsupported"},
    {4, 3, interaction_expansion_choice::unsupported,         "FLAVORS=4,SITES=3 -> unsupported"},
  };

  int fails = 0;
  for (const auto& c : cases) {
    const interaction_expansion_choice got =
        select_interaction_expansion(c.flavors, c.sites);
    if (got != c.expect) {
      std::printf("FAIL: %s (got %d, expected %d)\n", c.desc,
                  static_cast<int>(got), static_cast<int>(c.expect));
      ++fails;
    }
  }
  if (fails == 0) {
    std::printf("OK: all %zu interaction-expansion selection cases\n",
                sizeof(cases) / sizeof(cases[0]));
  }
  return fails ? 1 : 0;
}
