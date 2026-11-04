// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/input/looper.hpp>
#include <alps/input/parapack.hpp>
#include <alps/parameter_input.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

alps::Parameters text(std::string const & parameters) {
    if (parameters.empty())
        return alps::Parameters();
    std::istringstream in(parameters);
    return alps::read_parameter_text(in);
}

TEST(ParapackInput, StepsDefaultThermalizationToAnEighthOfTheSweeps) {
    alps::input::mc_steps_parameters const defaults = alps::input::read_mc_steps_parameters(text(""));
    EXPECT_EQ(defaults.sweeps.min(), 65536u);
    EXPECT_EQ(defaults.thermalization, 8192u);
    alps::input::mc_steps_parameters const p =
        alps::input::read_mc_steps_parameters(text("L = 100\nSWEEPS = 10*L\nTHERMALIZATION = SWEEPS/8"));
    EXPECT_EQ(p.sweeps.min(), 1000u);
    EXPECT_EQ(p.sweeps.max(), 1000u);
    EXPECT_EQ(p.thermalization, 125u);
    EXPECT_EQ(alps::input::read_mc_steps_parameters(text("SWEEPS = \"[100:200]\"")).sweeps.max(), 200u);
}

TEST(ParapackInput, AlgorithmFallsBackToWorker) {
    alps::input::parapack_algorithm_parameters const none = alps::input::read_parapack_algorithm_parameters(text(""));
    EXPECT_EQ(none.algorithm, "");
    EXPECT_FALSE(none.obsolete_worker_key);
    EXPECT_EQ(none.evaluator, "");
    alps::input::parapack_algorithm_parameters const worker =
        alps::input::read_parapack_algorithm_parameters(text("WORKER = \"ising\"\nEVALUATOR = \"default\""));
    EXPECT_EQ(worker.algorithm, "ising");
    EXPECT_TRUE(worker.obsolete_worker_key);
    EXPECT_EQ(worker.evaluator, "default");
    alps::input::parapack_algorithm_parameters const both =
        alps::input::read_parapack_algorithm_parameters(text("ALGORITHM = \"xy\"\nWORKER = \"ising\""));
    EXPECT_EQ(both.algorithm, "xy");
    EXPECT_FALSE(both.obsolete_worker_key);
}

TEST(ParapackInput, TaskAndSeeds) {
    alps::input::parapack_task_parameters const task = alps::input::read_parapack_task_parameters(text(""));
    EXPECT_EQ(task.num_clones.min(), 1u);
    EXPECT_EQ(task.num_clones.max(), 1u);
    EXPECT_FALSE(task.evaluate_only_finished_clones);
    EXPECT_TRUE(alps::input::read_parapack_task_parameters(text("EVALUATE_ONLY_FINISHED_CLONES = 1"))
        .evaluate_only_finished_clones);
    EXPECT_EQ(alps::input::read_parapack_task_parameters(text("N = 4\nNUM_CLONES = N")).num_clones.max(), 4u);

    EXPECT_FALSE(alps::input::read_parapack_clone_seeds(text("DISORDER_SEED = 5")).disorder_seed);
    alps::input::parapack_clone_seeds const seeds =
        alps::input::read_parapack_clone_seeds(text("SEED = 7\nDISORDER_SEED = 5"));
    EXPECT_EQ(*seeds.seed, 7u);
    EXPECT_EQ(*seeds.disorder_seed, 5u);

    alps::input::parapack_rng_parameters const rng =
        alps::input::read_parapack_rng_parameters(text("WORKER_SEED = 3\nDISORDER_SEED = 4"));
    EXPECT_EQ(rng.rng, "mt19937");
    EXPECT_EQ(rng.worker_seed, 3u);
    EXPECT_EQ(rng.disorder_seed, 4u);
    EXPECT_THROW(alps::input::read_parapack_rng_parameters(text("DISORDER_SEED = 4")), std::runtime_error);
}

TEST(ParapackInput, ExchangeTemperatures) {
    alps::input::exchange_temperature_parameters const set =
        alps::input::read_exchange_temperature_parameters(text("TEMPERATURE_SET = \"2 4\"\nBETA_MAX = 1\nBETA_MIN = 0.1"));
    EXPECT_EQ(set.distribution, 0);
    ASSERT_TRUE(set.inverse_temperatures);
    EXPECT_EQ(*set.inverse_temperatures, (std::vector<double>{0.5, 0.25}));

    alps::input::exchange_temperature_parameters const range =
        alps::input::read_exchange_temperature_parameters(text("NUM_REPLICAS = 4\nT_MAX = 2\nT_MIN = 0.5"));
    EXPECT_EQ(range.distribution, 2);
    EXPECT_FALSE(range.inverse_temperatures);
    EXPECT_EQ(range.num_replicas, 4);
    EXPECT_EQ(range.beta_max, 2.);
    EXPECT_EQ(range.beta_min, 0.5);

    EXPECT_THROW(alps::input::read_exchange_temperature_parameters(text("NUM_REPLICAS = 4")), std::invalid_argument);
    EXPECT_THROW(alps::input::read_exchange_temperature_parameters(
        text("TEMPERATURE_DISTRIBUTION_TYPE = 1\nNUM_REPLICAS = 4")), std::invalid_argument);
    EXPECT_THROW(alps::input::read_exchange_temperature_parameters(text("BETA_MAX = 1\nBETA_MIN = 0.1")),
                 std::invalid_argument);
}

TEST(ParapackInput, ExchangeOptimizationDefaultsFollowItsType) {
    alps::input::exchange_parameters const plain = alps::input::read_exchange_parameters(text(""));
    EXPECT_TRUE(plain.exchange);
    EXPECT_EQ(plain.interval, 1u);
    EXPECT_FALSE(plain.random_exchange);
    EXPECT_FALSE(plain.optimize);

    alps::input::exchange_parameters const rate =
        alps::input::read_exchange_parameters(text("OPTIMIZE_TEMPERATURE = 1\nSWEEPS = 800"));
    EXPECT_EQ(rate.optimization_type, alps::input::exchange_parameters::rate);
    EXPECT_EQ(rate.block_sweep_factor, 1.);
    EXPECT_EQ(rate.optimization_iterations, 1u);
    EXPECT_EQ(rate.initial_block_sweeps, 100);

    alps::input::exchange_parameters const population =
        alps::input::read_exchange_parameters(text("TEMPERATURE_OPTIMIZATION = 1\nOPTIMIZATION_TYPE = \"population\""));
    EXPECT_EQ(population.optimization_type, alps::input::exchange_parameters::population);
    EXPECT_EQ(population.block_sweep_factor, 2.);
    EXPECT_EQ(population.optimization_iterations, 7u);
    EXPECT_EQ(population.initial_block_sweeps, 512);

    EXPECT_THROW(alps::input::read_exchange_parameters(text("OPTIMIZE_TEMPERATURE = 1\nOPTIMIZATION_TYPE = \"other\"")),
                 std::runtime_error);
    EXPECT_THROW(alps::input::read_exchange_parameters(text("OPTIMIZE_TEMPERATURE = 1\nNUM_CLONES = 2")),
                 std::logic_error);
    // without exchange the optimization keys are not read
    EXPECT_NO_THROW(alps::input::read_exchange_parameters(
        text("NO_EXCHANGE = 1\nOPTIMIZE_TEMPERATURE = 1\nOPTIMIZATION_TYPE = \"other\"")));
}

TEST(ParapackInput, TemperatureScanAndWangLandau) {
    EXPECT_THROW(alps::input::read_temperature_scan_parameters(text("")), std::runtime_error);
    alps::input::temperature_scan_parameters const scan = alps::input::read_temperature_scan_parameters(
        text("NUM_TEMPERATURES = 3\nINITIAL_TEMPERATURE = 2\nDIFF_TEMPERATURE = -0.5"));
    EXPECT_EQ(scan.num_temperatures, 3);
    EXPECT_EQ(scan.initial_temperature, 2.);
    EXPECT_EQ(scan.diff_temperature, -0.5);
    EXPECT_FALSE(scan.initial_thermalization);

    alps::input::wanglandau_parameters const wl = alps::input::read_wanglandau_parameters(text(""));
    EXPECT_EQ(wl.check_interval, 16384u);
    EXPECT_EQ(wl.initial_update_factor, std::exp(1.0));
    EXPECT_EQ(wl.final_update_factor, 1.00000001);
    EXPECT_EQ(wl.flatness_threshold, 0.95);

    alps::Parameters const names = text("DIR_NAME = \".\"\nBASE_NAME = \"task1\"");
    EXPECT_EQ(alps::input::read_wanglandau_observable_dump(names).file, "task1.observable");
    EXPECT_EQ(alps::input::read_wanglandau_weight_dump(names).files, (std::vector<std::string>{"task1.weight"}));
    // each reads BASE_NAME only for its own default
    alps::input::wanglandau_weight_dump const weights =
        alps::input::read_wanglandau_weight_dump(text("DIR_NAME = \".\"\nWEIGHT_DUMP_FILE = \"a.w b.w\""));
    EXPECT_EQ(weights.files, (std::vector<std::string>{"a.w", "b.w"}));
}

TEST(LooperInput, TemperatureAndAnnealing) {
    EXPECT_FALSE(alps::input::read_looper_temperature_parameters(text("")).temperature);
    EXPECT_EQ(*alps::input::read_looper_temperature_parameters(text("BETA = 4")).temperature, 0.25);
    EXPECT_EQ(*alps::input::read_looper_temperature_parameters(text("T = 0.5\nBETA = 4")).temperature, 0.5);
    alps::input::looper_temperature_parameters const annealed = alps::input::read_looper_temperature_parameters(
        text("T = 0.1\nT_DURATION_0 = 100\nT_START_0 = 2\nBETA_DURATION_1 = 50\nBETA_START_1 = 2\nT_START_2 = 1"));
    EXPECT_EQ(*annealed.temperature, 0.1);
    ASSERT_EQ(annealed.annealing.size(), 2u);
    EXPECT_EQ(annealed.annealing[0], std::make_pair(100, 2.));
    EXPECT_EQ(annealed.annealing[1], std::make_pair(50, 0.5));
}

TEST(LooperInput, EstimatorsWeightsAndModel) {
    EXPECT_FALSE(alps::input::read_looper_worker_parameters(text("")).disable_improved_estimator);
    EXPECT_TRUE(alps::input::read_looper_worker_parameters(text("DISABLE_IMPROVED_ESTIMATOR = 0")).disable_improved_estimator);

    alps::input::looper_correlation_parameters const c =
        alps::input::read_looper_correlation_parameters(text("MEASURE[Green Function] = 1"));
    EXPECT_FALSE(c.measure_correlations);
    EXPECT_TRUE(c.measure_green_function);
    EXPECT_FALSE(c.measure_structure_factor);
    EXPECT_FALSE(alps::input::read_looper_initial_site(text("")));
    EXPECT_EQ(*alps::input::read_looper_initial_site(text("INITIAL_SITE = 3")), 3);

    alps::input::looper_weight_parameters const w = alps::input::read_looper_weight_parameters(text(""));
    EXPECT_FALSE(w.force_scatter);
    EXPECT_EQ(w.scattering_ratio, 0.5);
    EXPECT_EQ(*alps::input::read_looper_weight_parameters(text("FORCE_SCATTER = 0.2")).force_scatter, 0.2);

    alps::input::looper_model_parameters const m =
        alps::input::read_looper_model_parameters(text("USE_BOND_INDICES_AS_TYPES = 1\nLOOPER_DEBUG[MODEL OUTPUT] = \"cerr\""));
    EXPECT_EQ(m.disorder_seed, 0u);
    EXPECT_FALSE(m.use_site_indices_as_types);
    EXPECT_TRUE(m.use_bond_indices_as_types);
    EXPECT_EQ(*m.debug_model_output, "cerr");
}

}
