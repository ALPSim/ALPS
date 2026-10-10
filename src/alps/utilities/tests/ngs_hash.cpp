/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/hash.hpp>
#include <alps/ngs/stacktrace.hpp>

#include <cmath>
#include <stdio.h>
#include <iostream>
#include <stdexcept>

// TODO: check if matrix is linear independant (mathematica)
// tODO: improve hash
// after 2^32 change matrix
// if many independant generators are needed, take different matrecs. Differen meens diferent eigenvalues
class rng {
    public:
        rng(boost::uint64_t seed = 42)
            : state(seed)
        {
            if (seed == 0)
                throw std::runtime_error("Seed 0 is not valid" + ALPS_STACKTRACE);
        }
        boost::uint64_t operator()() {
            using alps::hash_value;
            return state = hash_value(state);
        }
    private:
        boost::uint64_t state;
};

#include <gtest/gtest.h>

TEST(HashGenerator, SeedAndSequenceAreReproducible) {
    EXPECT_THROW(rng(0), std::runtime_error);
    rng gen(42), repeated(42);
    for (std::size_t i = 0; i < 100000; ++i)
        ASSERT_EQ(gen(), repeated()) << "sequence index=" << i;
}

TEST(Hash, ZeroAndDistinctInputs) {
    EXPECT_EQ(alps::hash_value(boost::uint64_t(0)), 0u);
    EXPECT_NE(alps::hash_value(boost::uint64_t(1)), alps::hash_value(boost::uint64_t(2)));
}
