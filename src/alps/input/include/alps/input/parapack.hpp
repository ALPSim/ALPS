// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef ALPS_INPUT_PARAPACK_HPP
#define ALPS_INPUT_PARAPACK_HPP

#include <alps/parameter.h>
#include <alps/config.h>
#include <alps/parapack/integer_range.h>

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace alps { namespace input {

    // The Monte Carlo steps of a parapack worker.
    struct mc_steps_parameters {
        integer_range<unsigned int> sweeps;      // SWEEPS, default [65536:]
        unsigned int thermalization;             // THERMALIZATION, at least 0; default min(SWEEPS)/8
    };

    ALPS_DECL mc_steps_parameters read_mc_steps_parameters(Parameters const & parameters);

    // Which worker and evaluator parapack creates for a task.
    struct parapack_algorithm_parameters {
        std::string algorithm;                   // ALGORITHM, or the obsolete WORKER; empty if neither
        bool obsolete_worker_key;                // the name came from WORKER
        std::string evaluator;                   // EVALUATOR; empty if not set
    };

    ALPS_DECL parapack_algorithm_parameters read_parapack_algorithm_parameters(Parameters const & parameters);

    // NUM_CLONES, default 1: the range of clones a task runs.
    ALPS_DECL integer_range<std::uint32_t> read_parapack_num_clones(Parameters const & parameters);

    // EVALUATE_ONLY_FINISHED_CLONES, default false: leave running and
    // suspended clones out of the evaluation.
    ALPS_DECL bool read_evaluate_only_finished_clones(Parameters const & parameters);

    // The seeds a clone derives its worker and disorder seeds from.
    struct parapack_clone_seeds {
        std::optional<std::uint32_t> seed;          // SEED
        std::optional<std::uint32_t> disorder_seed; // DISORDER_SEED, read only with SEED
    };

    ALPS_DECL parapack_clone_seeds read_parapack_clone_seeds(Parameters const & parameters);

    // The random number generator of a worker. WORKER_SEED and DISORDER_SEED
    // are set by parapack for each clone.
    struct parapack_rng_parameters {
        std::string rng;                         // RNG, default mt19937
        std::uint32_t worker_seed;               // WORKER_SEED
        std::uint32_t disorder_seed;             // DISORDER_SEED
    };

    ALPS_DECL parapack_rng_parameters read_parapack_rng_parameters(Parameters const & parameters);

    // The inverse temperatures of the replicas in replica exchange.
    struct exchange_temperature_parameters {
        // 0: explicit, 1: equidistant in beta, 2: in 1/beta, 3: in sqrt(beta).
        // TEMPERATURE_DISTRIBUTION_TYPE; by default 1 with BETA_MAX and
        // BETA_MIN, 2 with T_MAX and T_MIN, and 0 with an explicit set.
        int distribution;
        // INVERSE_TEMPERATURE_SET, or the inverse of TEMPERATURE_SET
        std::optional<std::vector<double>> inverse_temperatures;
        // without an explicit set:
        int num_replicas;                        // NUM_REPLICAS
        double beta_max;                         // BETA_MAX, or 1/T_MIN
        double beta_min;                         // BETA_MIN, or 1/T_MAX
    };

    ALPS_DECL exchange_temperature_parameters read_exchange_temperature_parameters(Parameters const & parameters);

    // The replica exchange schedule.
    struct exchange_parameters {
        enum optimization { rate, population };

        bool exchange;                           // not NO_EXCHANGE, default true
        unsigned int interval;                   // EXCHANGE_INTERVAL, default 1
        bool random_exchange;                    // RANDOM_EXCHANGE, default false
        bool optimize;                           // OPTIMIZE_TEMPERATURE or TEMPERATURE_OPTIMIZATION, default false
        mc_steps_parameters steps;
        // read only when exchanging with temperature optimization:
        optimization optimization_type;          // OPTIMIZATION_TYPE: "rate" (default) or "population"
        double block_sweep_factor;               // BLOCK_SWEEP_FACTOR, default 1 (rate) or 2 (population)
        unsigned int optimization_iterations;    // OPTIMIZATION_ITERATIONS, default 1 (rate) or 7 (population)
        int initial_block_sweeps;                // INITIAL_BLOCK_SWEEPS, default THERMALIZATION (rate) or 512 (population)
    };

    ALPS_DECL exchange_parameters read_exchange_parameters(Parameters const & parameters);

    // PROCESS_PER_WORKER, default 1: processes that share one walker in
    // the parallel replica exchange worker.
    ALPS_DECL int read_process_per_worker(Parameters const & parameters);

    // The temperature scan of temperature_scan_adaptor.
    struct temperature_scan_parameters {
        int num_temperatures;                    // NUM_TEMPERATURES
        double initial_temperature;              // INITIAL_TEMPERATURE
        double diff_temperature;                 // DIFF_TEMPERATURE
        std::optional<int> initial_thermalization; // INITIAL_THERMALIZATION, at least 0: thermalization at the first temperature
    };

    ALPS_DECL temperature_scan_parameters read_temperature_scan_parameters(Parameters const & parameters);

    // The Wang-Landau learning stage.
    struct wanglandau_parameters {
        unsigned int check_interval;             // CHECK_INTERVAL, default 16384
        double initial_update_factor;            // INITIAL_UPDATE_FACTOR, default e
        double final_update_factor;              // FINAL_UPDATE_FACTOR, default 1.00000001
        double flatness_threshold;               // FLATNESS_THRESHOLD, default 0.95
    };

    ALPS_DECL wanglandau_parameters read_wanglandau_parameters(Parameters const & parameters);

    // Where the Wang-Landau observables and weights are written, relative
    // to DIR_NAME. DIR_NAME and BASE_NAME are set by parapack.
    struct wanglandau_observable_dump {
        std::string dir_name;                    // DIR_NAME
        std::string file;                        // OBSERVABLE_DUMP_FILE, default BASE_NAME.observable
    };

    ALPS_DECL wanglandau_observable_dump read_wanglandau_observable_dump(Parameters const & parameters);

    struct wanglandau_weight_dump {
        std::string dir_name;                    // DIR_NAME
        std::vector<std::string> files;          // WEIGHT_DUMP_FILE, default BASE_NAME.weight
    };

    ALPS_DECL wanglandau_weight_dump read_wanglandau_weight_dump(Parameters const & parameters);

} }

#endif
