// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef ALPS_INPUT_DMFT_HPP
#define ALPS_INPUT_DMFT_HPP

#include <alps/parameter.h>
#include <alps/config.h>

#include <optional>
#include <string>
#include <vector>

namespace alps { namespace input {

    // The parameters the dmft self-consistency program and its built-in
    // parts read. The program changes some parameters while it runs (H
    // during the first iteration, EPS_<f> and EPSSQ_<f> from the band
    // structure, BASENAME) and passes them all to the impurity solver, so
    // each part reads its parameters when it is set up, as before.

    // FLAVORS (default 2) and SITES (default 1).
    struct dmft_sizes {
        int flavors;
        int sites;
    };

    ALPS_DECL dmft_sizes read_dmft_sizes(Parameters const & parameters);

    // Which self-consistency loop and impurity solver dmft runs.
    struct dmft_run_parameters {
        bool cluster_loop;                       // CLUSTER_LOOP is set (not supported here)
        bool omega_loop;                         // OMEGA_LOOP, default false: loop in imaginary time
        std::string solver;                      // SOLVER, required unless CLUSTER_LOOP is set
        bool sc_write_delta;                     // SC_WRITE_DELTA, default false
        bool antiferromagnet;                    // ANTIFERROMAGNET, default false
        bool symmetrization;                     // SYMMETRIZATION, default false
        bool semicircle_hilbert;                 // SEMICIRCLE_HILBERT, default false
        dmft_sizes sizes;
    };

    ALPS_DECL dmft_run_parameters read_dmft_run_parameters(Parameters const & parameters);

    enum class dmft_loop {
        imaginary_time,         // G0(tau) is iterated
        dressed_imaginary_time, // G(tau) is iterated, for solvers that take the hybridization
        matsubara               // G0(i omega_n) is iterated
    };

    // The self-consistency loop. Fields a loop does not read are left at
    // their defaults.
    struct dmft_loop_parameters {
        int n_tau;                                     // N
        int n_matsubara;                               // NMATSUBARA (Matsubara loop)
        dmft_sizes sizes;                              // SITES is read by the Matsubara loop only
        double beta;                                   // BETA
        double h;                                      // H, default 0
        double mu;                                     // MU (not read by the dressed loop)
        double converged;                              // CONVERGED: largest change for convergence
        bool symmetrization;                           // SYMMETRIZATION
        double relax_rate;                             // RELAX_RATE, default 1 (Matsubara loop)
        int max_it;                                    // MAX_IT, default 1000
        std::string basename;                          // BASENAME
        std::optional<std::string> h_init;             // H_INIT: H for the initial G0
        std::string g0_tau_output;                     // G0TAU_output, default "G0tau_output"
        std::string g0_omega_output;                   // G0OMEGA_output, default "G0omega_output"
        std::optional<std::string> initial_g0_tau_output;  // G0TAU_input (Matsubara loop)
    };

    ALPS_DECL dmft_loop_parameters read_dmft_loop_parameters(Parameters const & parameters, dmft_loop loop);

    // The initial G0 of the Hilbert transformers.
    struct dmft_initial_g0_parameters {
        int n_tau;                                     // N
        int n_matsubara;                               // NMATSUBARA (Matsubara)
        int flavors;                                   // FLAVORS, default 2
        double beta;                                   // BETA (Matsubara)
        double mu;                                     // MU (Matsubara)
        double h;                                      // H, default 0 (Matsubara)
        bool insulating;                               // INSULATING is set: start from the atomic limit (Matsubara)
        std::optional<std::string> input;              // G0TAU_INPUT or G0OMEGA_INPUT, when not empty: read G0 from it
        std::optional<std::string> output;             // G0TAU_input or G0OMEGA_input: write the initial G0 to it
    };

    ALPS_DECL dmft_initial_g0_parameters read_dmft_initial_g0_tau(Parameters const & parameters);
    ALPS_DECL dmft_initial_g0_parameters read_dmft_initial_g0_omega(Parameters const & parameters);

    // ANTIFERROMAGNET, default false.
    ALPS_DECL bool read_dmft_antiferromagnet(Parameters const & parameters);

    enum class dmft_bandstructure { semicircle, square, hexagonal, dos };

    // DOSFILE set: a density of states from that file; TWODBS = hexagonal:
    // the hexagonal lattice; any other TWODBS: the square lattice; otherwise
    // the semicircular density of states.
    struct dmft_bandstructure_choice {
        dmft_bandstructure kind;
        std::string dosfile;
    };

    ALPS_DECL dmft_bandstructure_choice read_dmft_bandstructure(Parameters const & parameters);

    // The semicircular density of states: the hopping of each band
    // (flavors 2b and 2b+1), t<b> or t. Throws for an odd FLAVORS.
    ALPS_DECL std::vector<double> read_semicircle_hoppings(Parameters const & parameters);

    // The square and hexagonal lattices: t (default 1), tprime (default 0)
    // and the k-grid size L (default 128).
    struct dmft_lattice_hoppings {
        double t;
        double tprime;
        int L;
    };

    ALPS_DECL dmft_lattice_hoppings read_lattice_hoppings(Parameters const & parameters);

    // The first and second moments of the density of states of each flavor,
    // EPS_<f> and EPSSQ_<f>, when set.
    struct dmft_band_moments {
        std::vector<std::optional<double> > eps;
        std::vector<std::optional<double> > epssq;
    };

    ALPS_DECL dmft_band_moments read_band_moments(Parameters const & parameters, int flavors);

    // The Fourier transformers. MU and U are read only when asked for.
    struct dmft_fourier_parameters {
        dmft_sizes sizes;
        double beta;                                   // BETA
        double mu;                                     // MU
        double h;                                      // H, default 0
        double U;                                      // U
        std::vector<double> eps;                       // EPS_<f>, default 0
        std::vector<double> epssq;                     // EPSSQ_<f>, default 1
    };

    ALPS_DECL dmft_fourier_parameters read_dmft_fourier_parameters(Parameters const & parameters,
                                                                   bool with_mu, bool with_U);

    // The interaction matrix: from the file U_MATRIX, or from U (one or two
    // flavors) or U, J and U' (default U-2J).
    struct dmft_interaction_parameters {
        dmft_sizes sizes;
        std::optional<std::string> u_matrix;           // U_MATRIX
        double U;
        double J;
        double Uprime;
    };

    ALPS_DECL dmft_interaction_parameters read_dmft_interaction(Parameters const & parameters);

    // The external solver handoff. Fields a call does not read are left at
    // their defaults.
    struct dmft_external_solver_parameters {
        int n_tau;                                     // N
        int n_matsubara;                               // NMATSUBARA (Matsubara)
        dmft_sizes sizes;
        std::optional<std::string> tmpname;            // TMPNAME: base name of the input and output files
        bool sc_write_delta;                           // SC_WRITE_DELTA, default false: pass the hybridization
        bool antiferromagnet;                          // ANTIFERROMAGNET, default false (imaginary time, with SC_WRITE_DELTA)
        double beta;                                   // BETA (with SC_WRITE_DELTA)
        double mu;                                     // MU (Matsubara, with SC_WRITE_DELTA)
        double h;                                      // H, default 0 (Matsubara, with SC_WRITE_DELTA)
    };

    ALPS_DECL dmft_external_solver_parameters read_dmft_external_solver(Parameters const & parameters, dmft_loop loop);

    // What the hybridization handoff converts: MU, H (default 0) and
    // FLAVORS (required).
    struct dmft_hybridization_handoff {
        double mu;
        double h;
        int flavors;
    };

    ALPS_DECL dmft_hybridization_handoff read_dmft_hybridization_handoff(Parameters const & parameters);

    // MAX_TIME: the seconds a built-in impurity solver runs, default
    // default_time.
    ALPS_DECL unsigned int read_dmft_max_time(Parameters const & parameters, unsigned int default_time);

} }

#endif
