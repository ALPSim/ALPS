// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>
#include <alps/parapack/rng_helper.h>
#include <array>
#include <stdexcept>
#include <omp.h>

TEST(OpenMPWorkers, IndependentRandomEngines) {
    omp_set_dynamic(0);
    alps::Parameters parameters;
    parameters["WORKER_SEED"] = 42;
    parameters["DISORDER_SEED"] = 0;
    alps::rng_helper random(parameters);
    std::array<double, 2> samples{};
    int workers = 0;
#pragma omp parallel reduction(+:workers)
    {
        ++workers;
        if (omp_get_thread_num() < 2)
            samples[omp_get_thread_num()] = random.random_01(omp_get_thread_num());
    }
    ASSERT_EQ(workers, 2);
    EXPECT_NE(&random.engine(0), &random.engine(1));
    for (double sample : samples) {
        EXPECT_GE(sample, 0.0);
        EXPECT_LT(sample, 1.0);
    }
}
