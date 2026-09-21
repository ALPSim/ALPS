// Copyright (C) 2010-2012 by Lukas Gamper
//               2026      by the ALPS collaboration
// SPDX-License-Identifier: MIT

#include "ising.hpp"

#include <pyalps/export_simulation.hpp>
#include <nanobind/nanobind.h>

NB_MODULE(ising_c, m) {
    ALPS_EXPORT_SIM_TO_PYTHON(sim, ising_sim);
}
