// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef ALPS_INPUT_LOOPER_HPP
#define ALPS_INPUT_LOOPER_HPP

#include <alps/parameter.h>
#include <alps/config.h>

#include <optional>
#include <string>
#include <utility>
#include <vector>

// The parameters the loop program reads itself. The lattice and the model
// are read by the lattice and model libraries, and SWEEPS and THERMALIZATION
// through alps/input/parapack.hpp.
namespace alps { namespace input {

    // The temperature, and the annealing schedule that leads to it.
    struct looper_temperature_parameters {
        // T, TEMPERATURE, 1/BETA or 1/INVERSE_TEMPERATURE, in that order;
        // unset if none is given.
        std::optional<double> temperature;
        // Annealing stages i = 0, 1, ...: (duration in sweeps, start
        // temperature) from T_DURATION_i and T_START_i, or the TEMPERATURE_,
        // BETA_ or INVERSE_TEMPERATURE_ variants. The schedule ends at the
        // first stage without both keys, or with a duration or start that is
        // not positive. Read only with a temperature.
        std::vector<std::pair<int, double>> annealing;
    };

    ALPS_DECL looper_temperature_parameters read_looper_temperature_parameters(Parameters const & parameters);

    // What the SSE and path-integral workers read.
    struct looper_worker_parameters {
        bool disable_improved_estimator;         // DISABLE_IMPROVED_ESTIMATOR is set, whatever its value
    };

    ALPS_DECL looper_worker_parameters read_looper_worker_parameters(Parameters const & parameters);

    // Which correlation functions the correlation estimator measures.
    struct looper_correlation_parameters {
        bool measure_correlations;               // MEASURE[Correlations], default false
        bool measure_green_function;             // MEASURE[Green Function], default false
        bool measure_structure_factor;           // MEASURE[Structure Factor], default false
    };

    ALPS_DECL looper_correlation_parameters read_looper_correlation_parameters(Parameters const & parameters);

    // INITIAL_SITE: the site correlations are measured from; unset to
    // average over all sites.
    ALPS_DECL std::optional<int> read_looper_initial_site(Parameters const & parameters);

    // The vertex weights.
    struct looper_weight_parameters {
        // FORCE_SCATTER: the least fraction of scattering vertices. Unset,
        // it is 0.1 for frustrated models and 0 otherwise.
        std::optional<double> force_scatter;
        double scattering_ratio;                 // SCATTERING_RATIO, default 0.5
    };

    ALPS_DECL looper_weight_parameters read_looper_weight_parameters(Parameters const & parameters);

    // What loop reads about the model besides the model library.
    struct looper_model_parameters {
        unsigned int disorder_seed;              // DISORDER_SEED, default 0
        bool use_site_indices_as_types;          // USE_SITE_INDICES_AS_TYPES, default false
        bool use_bond_indices_as_types;          // USE_BOND_INDICES_AS_TYPES, default false
        // LOOPER_DEBUG[MODEL OUTPUT]: print the model, to standard error if
        // "cerr" and to standard output otherwise
        std::optional<std::string> debug_model_output;
    };

    ALPS_DECL looper_model_parameters read_looper_model_parameters(Parameters const & parameters);

    // The Wang-Landau steps of wl_steps.
    struct looper_wanglandau_parameters {
        bool multicanonical;                     // not DISABLE_MULTICANONICAL_MEASUREMENT
        bool zhou_bhatt;                         // not DISABLE_ZHOU_BHATT
        int block_sweeps;                        // BLOCK_SWEEPS, default 65536; without Zhou-Bhatt only
        int sweeps;                              // SWEEPS, default 65536
        int iterations;                          // WANG_LANDAU_ITERATIONS, default 16
    };

    ALPS_DECL looper_wanglandau_parameters read_looper_wanglandau_parameters(Parameters const & parameters);

} }

#endif
