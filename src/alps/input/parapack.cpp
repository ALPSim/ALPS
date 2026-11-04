// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/input/parapack.hpp>

#include <alps/expression.h>
#include <alps/numeric/double2int.hpp>
#include <alps/utilities/vectorio.hpp>
#include <boost/throw_exception.hpp>

#include <cmath>
#include <stdexcept>

namespace alps { namespace input {

    mc_steps_parameters read_mc_steps_parameters(Parameters const & p) {
        mc_steps_parameters result;
        result.sweeps = integer_range<unsigned int>(p.value_or_default("SWEEPS", "[65536:]"), p);
        result.thermalization = p.defined("THERMALIZATION")
            ? static_cast<int>(alps::evaluate("THERMALIZATION", p))
            : (result.sweeps.min BOOST_PREVENT_MACRO_SUBSTITUTION () >> 3);
        return result;
    }

    parapack_algorithm_parameters read_parapack_algorithm_parameters(Parameters const & p) {
        parapack_algorithm_parameters result;
        result.obsolete_worker_key = false;
        if (p.defined("ALGORITHM")) {
            result.algorithm = static_cast<std::string>(p["ALGORITHM"]);
        } else if (p.defined("WORKER")) {
            result.algorithm = static_cast<std::string>(p["WORKER"]);
            result.obsolete_worker_key = true;
        }
        if (p.defined("EVALUATOR"))
            result.evaluator = static_cast<std::string>(p["EVALUATOR"]);
        return result;
    }

    parapack_task_parameters read_parapack_task_parameters(Parameters const & p) {
        parapack_task_parameters result;
        result.num_clones = integer_range<std::uint32_t>(p.value_or_default("NUM_CLONES", "1"), p);
        result.evaluate_only_finished_clones = p.defined("EVALUATE_ONLY_FINISHED_CLONES") &&
            static_cast<bool>(alps::evaluate("EVALUATE_ONLY_FINISHED_CLONES", p));
        return result;
    }

    parapack_clone_seeds read_parapack_clone_seeds(Parameters const & p) {
        parapack_clone_seeds result;
        if (p.defined("SEED")) {
            result.seed = static_cast<std::uint32_t>(p["SEED"]);
            if (p.defined("DISORDER_SEED"))
                result.disorder_seed = static_cast<std::uint32_t>(p["DISORDER_SEED"]);
        }
        return result;
    }

    parapack_rng_parameters read_parapack_rng_parameters(Parameters const & p) {
        parapack_rng_parameters result;
        result.rng = static_cast<std::string>(p.value_or_default("RNG", "mt19937"));
        result.worker_seed = static_cast<std::uint32_t>(p["WORKER_SEED"]);
        result.disorder_seed = static_cast<std::uint32_t>(p["DISORDER_SEED"]);
        return result;
    }

    exchange_temperature_parameters read_exchange_temperature_parameters(Parameters const & p) {
        exchange_temperature_parameters result;
        result.distribution = -1;
        result.num_replicas = 0;
        result.beta_max = 0;
        result.beta_min = 0;
        if (p.defined("TEMPERATURE_DISTRIBUTION_TYPE")) {
            result.distribution = static_cast<int>(p["TEMPERATURE_DISTRIBUTION_TYPE"]);
        } else if (p.defined("BETA_MAX") && p.defined("BETA_MIN")) {
            result.distribution = 1;
        } else if (p.defined("T_MAX") && p.defined("T_MIN")) {
            result.distribution = 2;
        }
        if (p.defined("INVERSE_TEMPERATURE_SET") || p.defined("TEMPERATURE_SET")) result.distribution = 0;
        if (result.distribution < 0 || result.distribution > 3)
            boost::throw_exception(std::invalid_argument("illegal temperature distribution"));
        if (p.defined("INVERSE_TEMPERATURE_SET")) {
            std::vector<double> beta;
            read_vector_resize(p["INVERSE_TEMPERATURE_SET"], beta);
            result.inverse_temperatures = beta;
        } else if (p.defined("TEMPERATURE_SET")) {
            std::vector<double> beta;
            read_vector_resize(p["TEMPERATURE_SET"], beta);
            for (std::size_t i = 0; i < beta.size(); ++i) beta[i] = 1 / beta[i];
            result.inverse_temperatures = beta;
        } else if (p.defined("NUM_REPLICAS")) {
            result.num_replicas = static_cast<int>(evaluate("NUM_REPLICAS", p));
            if (p.defined("BETA_MAX") && p.defined("BETA_MIN")) {
                result.beta_max = evaluate("BETA_MAX", p);
                result.beta_min = evaluate("BETA_MIN", p);
            } else if (p.defined("T_MAX") && p.defined("T_MIN")) {
                result.beta_max = 1 / evaluate("T_MIN", p);
                result.beta_min = 1 / evaluate("T_MAX", p);
            } else {
                boost::throw_exception(std::invalid_argument("MAX/MIN values of BETA (or T) not defined"));
            }
        } else {
            boost::throw_exception(std::invalid_argument("temperature set not defined"));
        }
        return result;
    }

    exchange_parameters read_exchange_parameters(Parameters const & p) {
        exchange_parameters result;
        result.exchange = !static_cast<bool>(evaluate(p.value_or_default("NO_EXCHANGE", 0)));
        result.interval = static_cast<int>(evaluate(p.value_or_default("EXCHANGE_INTERVAL", 1)));
        result.random_exchange = static_cast<bool>(evaluate(p.value_or_default("RANDOM_EXCHANGE", 0)));
        result.optimize = static_cast<bool>(p.value_or_default("OPTIMIZE_TEMPERATURE", 0)) ||
                          static_cast<bool>(p.value_or_default("TEMPERATURE_OPTIMIZATION", 0));
        result.steps = read_mc_steps_parameters(p);
        result.optimization_type = exchange_parameters::rate;
        result.block_sweep_factor = 1;
        result.optimization_iterations = 0;
        result.initial_block_sweeps = 0;
        if (result.exchange && result.optimize) {
            if (!p.defined("OPTIMIZATION_TYPE") || p["OPTIMIZATION_TYPE"] == "rate") {
                result.optimization_type = exchange_parameters::rate;
            } else if (p["OPTIMIZATION_TYPE"] == "population") {
                result.optimization_type = exchange_parameters::population;
            } else {
                boost::throw_exception(
                    std::runtime_error("unknown OPIMIZATION_TYPE: " + p["OPTIMIZATION_TYPE"]));
            }
            if (p.defined("NUM_CLONES") &&
                integer_range<std::uint32_t>(p["NUM_CLONES"]).max BOOST_PREVENT_MACRO_SUBSTITUTION () > 1)
                boost::throw_exception(std::logic_error("temperature optimization and multipile clones"
                                                        " are exclusive"));
            if (result.optimization_type == exchange_parameters::rate) {
                result.block_sweep_factor = p.value_or_default("BLOCK_SWEEP_FACTOR", 1);
                result.optimization_iterations =
                    static_cast<unsigned int>(p.value_or_default("OPTIMIZATION_ITERATIONS", 1));
                result.initial_block_sweeps = p.value_or_default("INITIAL_BLOCK_SWEEPS", result.steps.thermalization);
            } else {
                result.block_sweep_factor = p.value_or_default("BLOCK_SWEEP_FACTOR", 2);
                result.optimization_iterations =
                    static_cast<unsigned int>(p.value_or_default("OPTIMIZATION_ITERATIONS", 7));
                result.initial_block_sweeps = p.value_or_default("INITIAL_BLOCK_SWEEPS", (2 << 8));
            }
        }
        return result;
    }

    int read_process_per_worker(Parameters const & p) {
        return p.value_or_default("PROCESS_PER_WORKER", 1);
    }

    temperature_scan_parameters read_temperature_scan_parameters(Parameters const & p) {
        temperature_scan_parameters result;
        result.num_temperatures = static_cast<int>(evaluate("NUM_TEMPERATURES", p));
        result.initial_temperature = evaluate("INITIAL_TEMPERATURE", p);
        result.diff_temperature = evaluate("DIFF_TEMPERATURE", p);
        if (p.defined("INITIAL_THERMALIZATION"))
            result.initial_thermalization = static_cast<int>(evaluate("INITIAL_THERMALIZATION", p));
        return result;
    }

    wanglandau_parameters read_wanglandau_parameters(Parameters const & p) {
        wanglandau_parameters result;
        result.check_interval = 16384;
        result.initial_update_factor = std::exp(1.0);
        result.final_update_factor = 1.00000001;
        result.flatness_threshold = 0.95;
        if (p.defined("CHECK_INTERVAL"))
            result.check_interval = numeric::double2int(evaluate("CHECK_INTERVAL", p));
        if (p.defined("INITIAL_UPDATE_FACTOR"))
            result.initial_update_factor = evaluate("INITIAL_UPDATE_FACTOR", p);
        if (p.defined("FINAL_UPDATE_FACTOR"))
            result.final_update_factor = evaluate("FINAL_UPDATE_FACTOR", p);
        if (p.defined("FLATNESS_THRESHOLD"))
            result.flatness_threshold = evaluate("FLATNESS_THRESHOLD", p);
        return result;
    }

    wanglandau_observable_dump read_wanglandau_observable_dump(Parameters const & p) {
        wanglandau_observable_dump result;
        result.dir_name = static_cast<std::string>(p["DIR_NAME"]);
        if (p.defined("OBSERVABLE_DUMP_FILE"))
            result.file = static_cast<std::string>(p["OBSERVABLE_DUMP_FILE"]);
        else
            result.file = static_cast<std::string>(p["BASE_NAME"]) + ".observable";
        return result;
    }

    wanglandau_weight_dump read_wanglandau_weight_dump(Parameters const & p) {
        wanglandau_weight_dump result;
        result.dir_name = static_cast<std::string>(p["DIR_NAME"]);
        if (p.defined("WEIGHT_DUMP_FILE"))
            read_vector_resize(p["WEIGHT_DUMP_FILE"], result.files);
        else
            result.files.push_back(p["BASE_NAME"] + ".weight");
        return result;
    }

} }
