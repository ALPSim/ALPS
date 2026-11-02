// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef ALPS_INPUT_HYBRIDIZATION_HPP
#define ALPS_INPUT_HYBRIDIZATION_HPP

#include <alps/params.hpp>
#include <alps/config.h>

#include <string>

namespace alps { namespace input {

    // Defines every parameter the hybridization-expansion (CT-HYB) solver
    // reads, with its type and default. SEED is defined by alps::mcbase.
    ALPS_DECL void define_hybridization_parameters(params & parameters);

    // The values of the parameters define_hybridization_parameters defined.
    struct hybridization_parameters {
        // simulation
        long sweeps;                       // SWEEPS after thermalization
        long thermalization;               // THERMALIZATION sweeps
        long max_time;                     // MAX_TIME in seconds
        int n_meas;                        // N_MEAS: updates between measurements
        int output_period;                 // OUTPUT_PERIOD, default 100000 sweeps
        int n_histogram_orders;            // N_HISTOGRAM_ORDERS, default 50
        bool timeseries;                   // TIMESERIES, default false: fixed instead of detailed binning
        int num_bins;                      // NUM_BINS, default 0
        bool verbose;                      // VERBOSE, default false
        // updates
        bool spinflip;                     // SPINFLIP, default false
        bool globalflip;                   // GLOBALFLIP, default false
        // measurements
        bool measure_time;                 // MEASURE_time, default true
        bool measure_freq;                 // MEASURE_freq, default false
        bool measure_legendre;             // MEASURE_legendre, default false
        bool measure_nn;                   // MEASURE_nn, default false
        bool measure_nnt;                  // MEASURE_nnt, default false
        bool measure_nnw;                  // MEASURE_nnw, default false
        bool measure_g2w;                  // MEASURE_g2w, default false
        bool measure_h2w;                  // MEASURE_h2w, default false
        bool measure_sector_statistics;    // MEASURE_sector_statistics, default false
        bool compute_vertex;               // COMPUTE_VERTEX, default false
        bool accurate_covariance;          // ACCURATE_COVARIANCE, default true
        int n_matsubara;                   // N_MATSUBARA, default 0
        int n_legendre;                    // N_LEGENDRE, default 0
        int n_nn;                          // N_nn, default 0
        int n_W;                           // N_W: bosonic frequencies, default 0
        int n_w2;                          // N_w2: fermionic frequencies of two-particle functions, default 0
        // input and output
        std::string basepath;              // BASEPATH, default ""
        bool text_output;                  // TEXT_OUTPUT, default false
        bool delta_in_hdf5;                // DELTA_IN_HDF5, default false
        bool dmft_framework;               // DMFT_FRAMEWORK, default false
        std::string u_matrix;              // U_MATRIX, default "": U, U' and J instead
        bool umatrix_in_hdf5;              // UMATRIX_IN_HDF5, default false
        std::string mu_vector;             // MU_VECTOR, default "": MU instead
        bool mu_in_hdf5;                   // MU_IN_HDF5, default false
        std::string ret_int_k;             // RET_INT_K, default "": no retarded interaction
        bool k_in_hdf5;                    // K_IN_HDF5, default false
        // physics
        double beta;                       // BETA
        int n_tau;                         // N_TAU
        int n_orbitals;                    // N_ORBITALS
        std::string delta;                 // DELTA: file with the hybridization function
        double mu;                         // MU, default 0
        double J;                          // J, default 0
        double U;                          // U; 0 when U_MATRIX is given
        double U_prime;                    // U', default U-2J; 0 when U_MATRIX is given
        bool u_ignored;                    // U is given, and U_MATRIX replaces it
        bool mu_ignored;                   // MU is given, and MU_VECTOR replaces it
    };

    // Reads parameters on which define_hybridization_parameters has been
    // called, and checks the ones that depend on each other.
    // @throws std::invalid_argument for a measurement without its sizes
    ALPS_DECL hybridization_parameters read_hybridization_parameters(params const & parameters);

} }

#endif
