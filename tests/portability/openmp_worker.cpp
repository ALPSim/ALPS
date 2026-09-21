// SPDX-License-Identifier: MIT
#include <alps/parapack/rng_helper.h>
#include <array>
#include <stdexcept>
#include <omp.h>

int main() {
    alps::Parameters parameters;
    parameters["WORKER_SEED"] = 42;
    parameters["DISORDER_SEED"] = 0;
    alps::rng_helper random(parameters);
    std::array<double, 2> samples{};
    int workers = 0;
#pragma omp parallel reduction(+:workers)
    {
        ++workers;
        samples[omp_get_thread_num()] = random.random_01(omp_get_thread_num());
    }
    if (workers != 2 || &random.engine(0) == &random.engine(1))
        throw std::runtime_error("OpenMP workers need independent random engines");
    for (double sample : samples)
        if (!(sample >= 0.0 && sample < 1.0))
            throw std::runtime_error("Invalid worker random sample");
}
