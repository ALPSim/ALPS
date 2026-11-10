// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef ALPS_INPUT_HYBRIDIZATION_HPP
#define ALPS_INPUT_HYBRIDIZATION_HPP

#include <alps/params.hpp>
#include <alps/config.h>

#include <cstddef>
#include <cstdint>
#include <string>

namespace alps { namespace input {

    // Defines every parameter the hybridization-expansion (CT-HYB) solver
    // reads, with its type and default. SEED is defined by alps::mcbase.
    ALPS_DECL void define_hybridization_parameters(params & parameters);

    // The values of the parameters define_hybridization_parameters defined.
    // Sweep counts and sizes are unsigned, so a negative value is rejected
    // when it is read.
    struct hybridization_parameters {
        // simulation
        std::uint64_t sweeps;              // SWEEPS after thermalization
        std::uint64_t thermalization;      // THERMALIZATION sweeps
        long max_time;                     // MAX_TIME in seconds
        std::size_t n_meas;                // N_MEAS: updates between measurements
        int output_period;                 // OUTPUT_PERIOD
        std::size_t n_histogram_orders;    // N_HISTOGRAM_ORDERS
        bool timeseries;                   // TIMESERIES: fixed instead of detailed binning
        std::size_t num_bins;              // NUM_BINS
        bool verbose;                      // VERBOSE
        // updates
        bool spinflip;                     // SPINFLIP
        bool globalflip;                   // GLOBALFLIP
        // measurements
        bool measure_time;                 // MEASURE_time
        bool measure_freq;                 // MEASURE_freq
        bool measure_legendre;             // MEASURE_legendre
        bool measure_nn;                   // MEASURE_nn
        bool measure_nnt;                  // MEASURE_nnt
        bool measure_nnw;                  // MEASURE_nnw
        bool measure_g2w;                  // MEASURE_g2w
        bool measure_h2w;                  // MEASURE_h2w
        bool measure_sector_statistics;    // MEASURE_sector_statistics
        bool compute_vertex;               // COMPUTE_VERTEX
        bool accurate_covariance;          // ACCURATE_COVARIANCE
        std::size_t n_matsubara;           // N_MATSUBARA
        std::size_t n_legendre;            // N_LEGENDRE
        std::size_t n_nn;                  // N_nn
        std::size_t n_W;                   // N_W: bosonic frequencies
        std::size_t n_w2;                  // N_w2: fermionic frequencies of two-particle functions
        // input and output
        std::string basepath;              // BASEPATH
        bool text_output;                  // TEXT_OUTPUT
        bool delta_in_hdf5;                // DELTA_IN_HDF5
        bool dmft_framework;               // DMFT_FRAMEWORK
        std::string u_matrix;              // U_MATRIX; empty: U, U' and J instead
        bool umatrix_in_hdf5;              // UMATRIX_IN_HDF5
        std::string mu_vector;             // MU_VECTOR; empty: MU instead
        bool mu_in_hdf5;                   // MU_IN_HDF5
        std::string ret_int_k;             // RET_INT_K; empty: no retarded interaction
        bool k_in_hdf5;                    // K_IN_HDF5
        // physics
        double beta;                       // BETA
        std::size_t n_tau;                 // N_TAU
        std::size_t n_orbitals;            // N_ORBITALS
        std::string delta;                 // DELTA: file with the hybridization function
        double mu;                         // MU
        double J;                          // J
        double U;                          // U; 0 when U_MATRIX is given
        double U_prime;                    // U'; 0 when U_MATRIX is given
        bool u_ignored;                    // U is given, and U_MATRIX replaces it
        bool mu_ignored;                   // MU is given, and MU_VECTOR replaces it
    };

    // Reads parameters on which define_hybridization_parameters has been
    // called, and checks the ones that depend on each other.
    // @throws std::invalid_argument for a measurement without its sizes
    // @throws alps::params_ns::exception::value_mismatch for a negative count or size
    ALPS_DECL hybridization_parameters read_hybridization_parameters(params const & parameters);

} }

#endif
