// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// The random streams of the MPI ranks of one job do not repeat those of a job
// with a neighbouring SEED. Rank k used SEED + k, so rank 1 of SEED=42 drew
// the numbers of rank 0 of SEED=43, and averaging such jobs under-reported
// the error. Rank 0 keeps SEED itself.

#include <alps/mc/mcbase.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <gtest/gtest.h>

#include <vector>

namespace {
class draws : public alps::mcbase {
public:
    draws(parameters_type const& params, std::size_t rank) : alps::mcbase(params, rank) {}
    void update() {}
    void measure() {}
    double fraction_completed() const { return 1.; }
    std::vector<double> take(int n) {
        std::vector<double> x(n);
        for (double& v : x) v = random();
        return x;
    }
};

std::vector<double> stream(long seed, std::size_t rank) {
    draws::parameters_type params;
    draws::define_parameters(params);
    params["SEED"] = seed;
    return draws(params, rank).take(100);
}
}

TEST(mc, RankStreamsDoNotRepeatANeighbouringSeed) {
    // rank 0 is the plain Mersenne twister with SEED, as before
    boost::mt19937 engine(42);
    boost::uniform_01<double> uniform;
    std::vector<double> const rank0 = stream(42, 0);
    for (double x : rank0) EXPECT_EQ(x, uniform(engine));

    std::vector<double> const rank1 = stream(42, 1);
    EXPECT_NE(rank1, stream(43, 0));
    EXPECT_NE(stream(42, 2), stream(43, 1));
}

