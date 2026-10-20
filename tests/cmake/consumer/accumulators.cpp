// SPDX-License-Identifier: MIT
#include <alps/accumulators.hpp>
#include <alps/hdf5/archive.hpp>
#include <cstdio>
#include <stdexcept>

int main() {
    namespace acc = alps::accumulators;
    acc::accumulator_set measurements;
    measurements << acc::FullBinningAccumulator<double>("x");
    for (int i = 0; i < 1024; ++i)
        measurements["x"] << 0.5 * (i % 4);
    const acc::result_set results(measurements);
    {
        alps::hdf5::archive archive("accumulators_contract.h5", "w");
        archive["results"] << results;
    }
    acc::result_set loaded;
    {
        alps::hdf5::archive archive("accumulators_contract.h5", "r");
        archive["results"] >> loaded;
    }
    std::remove("accumulators_contract.h5");
    if (loaded["x"].count() != 1024 || loaded["x"].mean<double>() != results["x"].mean<double>()
        || loaded["x"].error<double>() != results["x"].error<double>())
        throw std::runtime_error("accumulators archive round trip failed");
}
