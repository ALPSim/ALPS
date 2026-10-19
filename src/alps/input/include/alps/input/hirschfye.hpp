// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef ALPS_INPUT_HIRSCHFYE_HPP
#define ALPS_INPUT_HIRSCHFYE_HPP

#include <alps/parameter.h>
#include <alps/config.h>

#include <optional>
#include <string>

namespace alps { namespace input {

    // The parameters the Hirsch-Fye impurity solver reads, built into dmft
    // or as the hirschfye program.
    struct hirschfye_parameters {
        int thermalization;                      // THERMALIZATION: sweeps before measuring
        int sweeps;                              // SWEEPS: measured sweeps
        double beta;                             // BETA
        int n_tau;                               // N: time slices
        int sites;                               // SITES, default 1
        double U;                                // U
        int fourpoint_interval;                  // FOURPOINT_INTERVAL, default 1000
        bool measure_fourpoint;                  // MEASURE_FOURPOINT_FUNCTION, default false
        double tolerance;                        // TOLERANCE: allowed loss of precision
        int max_time;                            // MAX_TIME in seconds, default 86400
        int n_matsubara;                         // NMATSUBARA
    };

    ALPS_DECL hirschfye_parameters read_hirschfye_parameters(Parameters const & parameters);

    // What the solver reads when it returns its result.
    struct hirschfye_result_parameters {
        int n_tau;                               // N
        int n_matsubara;                         // NMATSUBARA
        int sites;                               // SITES, required here
        std::optional<std::string> checkpoint;   // CHECKPOINT: write a checkpoint to this file
        std::optional<std::string> outfile;      // OUTFILE: the hirschfye program's output file
    };

    ALPS_DECL hirschfye_result_parameters read_hirschfye_result_parameters(Parameters const & parameters);

    // The hirschfye program's input: NMATSUBARA, SITES (default 1) and
    // FLAVORS (default 2).
    struct hirschfye_input_sizes {
        int n_matsubara;
        int sites;
        int flavors;
    };

    ALPS_DECL hirschfye_input_sizes read_hirschfye_input_sizes(Parameters const & parameters);

} }

#endif
