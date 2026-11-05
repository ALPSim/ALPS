// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef ALPS_INPUT_CTINT_HPP
#define ALPS_INPUT_CTINT_HPP

#include <alps/params.hpp>
#include <alps/config.h>

#include <cstdint>
#include <string>

namespace alps { namespace input {

    // Defines every parameter the interaction-expansion (CT-INT) solver
    // reads, with its type and default. SEED is defined by alps::mcbase.
    ALPS_DECL void define_ctint_parameters(params & parameters);

    // The values of the parameters define_ctint_parameters defined.
    struct ctint_parameters {
        unsigned int max_order;                  // MAX_ORDER
        unsigned int flavors;                    // FLAVORS, else N_ORBITALS
        unsigned int sites;                      // SITES
        unsigned int n_matsubara;                // N_MATSUBARA, or its older name NMATSUBARA
        unsigned int n_matsubara_measurements;   // NMATSUBARA_MEASUREMENTS, else n_matsubara
        unsigned int n_tau;                      // N_TAU, or its older name N
        unsigned int n_self;                     // NSELF, else 10 n_tau
        std::uint64_t sweeps;                    // SWEEPS: measurement steps
        unsigned long thermalization;            // THERMALIZATION: steps before measuring
        long max_time;                           // MAX_TIME in seconds
        double beta;                             // BETA
        double U;                                // U: on-site interaction
        double alpha;                            // ALPHA: shift of the auxiliary field
        unsigned int recalc_period;              // RECALC_PERIOD
        unsigned int measurement_period;         // MEASUREMENT_PERIOD, else N_MEAS
        unsigned int convergence_check_period;   // CONVERGENCE_CHECK_PERIOD, else RECALC_PERIOD
        double almost_zero;                      // ALMOSTZERO
        bool histogram_measurement;              // HISTOGRAM_MEASUREMENT
        bool atomic;                             // ATOMIC: start from the atomic limit
        bool dmft_framework;                     // DMFT_FRAMEWORK: G0 in INFILE as dmft writes it
        unsigned int tau_discretization_for_exp; // TAU_DISCRETIZATION_FOR_EXP; 0: exp(i omega tau) on the fly
        std::string infile;                      // INFILE with G0; empty if ATOMIC and not given
    };

    // Reads parameters on which define_ctint_parameters has been called.
    ALPS_DECL ctint_parameters read_ctint_parameters(params const & parameters);

} }

#endif
