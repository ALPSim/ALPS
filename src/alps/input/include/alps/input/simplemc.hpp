// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef ALPS_INPUT_SIMPLEMC_HPP
#define ALPS_INPUT_SIMPLEMC_HPP

#include <alps/parameter.h>
#include <alps/config.h>

#include <string>
#include <vector>

namespace alps { namespace input {

    // The parameters the Ising, XY and Heisenberg workers of simplemc read.
    // SWEEPS and THERMALIZATION are read through alps/input/parapack.hpp.
    struct simplemc_parameters {
        double beta;                             // 1/T, or 0 without T
        double field;                            // H, default 0
        std::vector<double> couplings;           // per bond type t: Jt, or J (default 1)
        int snapshot_interval;                   // SNAPSHOT_INTERVAL in sweeps, default 0: no snapshots
        std::string base_name;                   // BASE_NAME, set by parapack
        int clone_id;                            // CLONE_ID, set by parapack
    };

    ALPS_DECL simplemc_parameters read_simplemc_parameters(Parameters const & parameters, int num_bond_types);

} }

#endif
