/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef ALPS_UTILITIES_STREAM_SEED_HPP
#define ALPS_UTILITIES_STREAM_SEED_HPP

#include <cstdint>

namespace alps {

    /// Seed for random stream `stream` (an MPI rank or a run) of a simulation
    /// with seed `seed`.
    ///
    /// Stream 0 uses `seed` itself, so a single run is unchanged. The other
    /// streams get a scrambled seed in [0, 2^31) instead of `seed + stream`:
    /// with `seed + stream`, stream 1 of a job with SEED=42 was stream 0 of the
    /// job with SEED=43, so jobs with consecutive seeds shared random numbers.
    inline std::uint64_t stream_seed(std::uint64_t seed, std::uint64_t stream) {
        if (stream == 0)
            return seed;
        // splitmix64 of the seed, mixed with the stream and scrambled again
        auto mix = [](std::uint64_t z) {
            z += 0x9e3779b97f4a7c15ULL;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            return z ^ (z >> 31);
        };
        return mix(mix(seed) ^ stream) & 0x7fffffffULL;
    }

}

#endif
