/*****************************************************************************
 *
 * Copyright (C) 2005-2026 by Emanuel Gull
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include "selfconsistency.h"

void define_selfconsistency_parameters(alps::params &p){

  p.define<int>("N", "Number of imaginary time discretization points");
  p.define<int>("NMATSUBARA", "Number of matsubara frequency discretization points");
  p.define<int>("FLAVORS", 2, "Number of spins or (diagonal) orbitals");
  p.define<double>("BETA", "Inverse temperature");
  p.define<double>("MU", "Chemical Potential");
  p.define<double>("t", "Nearest-neighbor hopping");
  p.define<double>("tprime", 0., "Next-nearest-neighbor hopping");

  p.define<double>("sc.SELFCONSISTENCY_CONVERGENCE", 1.e-10, "Difference in Green's function at which the iteration is converged");
  p.define<int>("sc.MAX_IT", 1000, "Maximum number of DCA iterations");
  p.define<bool>("sc.ADJUST_MU", false, "set to true to find target density, rather than run for fixed mu");
  p.define<double>("sc.TARGET_DENSITY", "desired target density, see also sc.ADJUST_MU");
  p.define<bool>("sc.SPIN_SYMMETRIZATION", true, "symmetrize 'up' and 'down' spins after each solver call.");
}
