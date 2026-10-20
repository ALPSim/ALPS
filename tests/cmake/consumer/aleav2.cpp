// SPDX-License-Identifier: MIT
#include <alps/alea.hpp>
#include <alps/hdf5/archive.hpp>
#include <cstdio>
#include <stdexcept>

int main() {
    alps::alea::var_acc<double> accumulator(1);
    for (int i = 0; i < 100; ++i)
        accumulator << 0.5 * i;
    const alps::alea::var_result<double> result = accumulator.finalize();
    {
        alps::hdf5::archive archive("aleav2_contract.h5", "w");
        alps::alea::hdf5_serializer serializer(archive, "");
        alps::alea::serialize(serializer, "x", result);
    }
    alps::alea::var_result<double> loaded;
    {
        alps::hdf5::archive archive("aleav2_contract.h5", "r");
        alps::alea::hdf5_serializer serializer(archive, "");
        alps::alea::deserialize(serializer, "x", loaded);
    }
    std::remove("aleav2_contract.h5");
    if (loaded.count() != 100 || loaded.mean()[0] != result.mean()[0]
        || loaded.stderror()[0] != result.stderror()[0])
        throw std::runtime_error("ALEA archive round trip failed");
}
