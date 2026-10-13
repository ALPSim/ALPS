/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2003 by Brigitte Surer and Jan Gukelberger                        *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <alps/scheduler/montecarlo.h>
#include <alps/alea.h>
#include <boost/random.hpp>
#include <boost/multi_array.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>

class Simulation
{
public:
    Simulation(double beta,size_t L, std::string output_file)
    :   eng_(42)
    ,   rng_(eng_, dist_)
    ,   L_(L)
    ,   beta_(beta)
    ,   spins_(boost::extents[L][L])
    ,   energy_("E")
    ,   magnetization_("m")
    ,   abs_magnetization_("|m|")
    ,   m2_("m^2")
    ,   m4_("m^4")
    ,   filename_(output_file)
    {
        // Init exponential map
        for(int E = -4; E <= 4; E += 2)
            exp_table_[E] = exp(2*beta*E);

        // Init random spin configuration
        for(size_t i = 0; i < L; ++i)
        {
            for(size_t j = 0; j < L; ++j)
                spins_[i][j] = 2 * randint(2) - 1;
        }
    }

    void run(size_t ntherm,size_t n)
    {
        thermalization_ = ntherm;
        sweeps_=n;
        // Thermalize for ntherm steps
        while(ntherm--)
            step();

        // Run n steps
        while(n--)
        {
            step();
            measure();
        }

        //save the observables
        save(filename_);

    }

    void verify_archive() const {
        alps::hdf5::archive archive(filename_, "r");
        std::size_t lattice = 0, sweeps = 0, thermalization = 0;
        double beta = -1.;
        archive["/parameters/L"] >> lattice;
        archive["/parameters/BETA"] >> beta;
        archive["/parameters/SWEEPS"] >> sweeps;
        archive["/parameters/THERMALIZATION"] >> thermalization;
        EXPECT_EQ(lattice, L_);
        EXPECT_EQ(beta, beta_);
        EXPECT_EQ(sweeps, sweeps_);
        EXPECT_EQ(thermalization, thermalization_);
        for (const auto* observable : {&energy_, &magnetization_, &abs_magnetization_, &m2_, &m4_}) {
            SCOPED_TRACE(observable->name());
            const auto path = "/simulation/results/" + observable->representation();
            unsigned long long count = 0;
            double mean = 0., error = 0.;
            archive[path + "/count"] >> count;
            archive[path + "/mean/value"] >> mean;
            archive[path + "/mean/error"] >> error;
            EXPECT_EQ(count, sweeps_);
            EXPECT_EQ(mean, observable->mean());
            EXPECT_EQ(error, observable->error());
            EXPECT_TRUE(std::isfinite(mean));
            EXPECT_TRUE(std::isfinite(error));
            EXPECT_GE(error, 0.);
        }
        // Exact physical bounds; no probability of rejecting a valid sample.
        EXPECT_GE(energy_.mean(), -2.);
        EXPECT_LE(energy_.mean(), 2.);
        EXPECT_GE(magnetization_.mean(), -1.);
        EXPECT_LE(magnetization_.mean(), 1.);
        EXPECT_GE(abs_magnetization_.mean(), 0.);
        EXPECT_LE(abs_magnetization_.mean(), 1.);
        EXPECT_LE(m4_.mean(), m2_.mean());
        EXPECT_LE(m2_.mean(), abs_magnetization_.mean());
    }
    void step()
    {
        for(size_t s = 0; s < L_*L_; ++s)
        {
            // Pick random site k=(i,j)
            int i = randint(L_);
            int j = randint(L_);

            // Measure local energy e = -s_k * sum_{l nn k} s_l
            int e = spins_[(i-1+L_)%L_][j] + spins_[(i+1)%L_][j] +
            spins_[i][(j-1+L_)%L_] + spins_[i][(j+1)%L_];
            e *= -spins_[i][j];

            // Flip s_k with probability exp(2 beta e)
            if(e > 1 || rng_() < exp_table_[e])
                spins_[i][j] = -spins_[i][j];
        }
    }
    void measure()
    {
        int E = 0; // energy
        int M = 0; // magnetization
        for(size_t i = 0; i < L_; ++i)
        {
            for(size_t j = 0; j < L_; ++j)
            {
                E -= spins_[i][j]*(spins_[(i+1)%L_][j] + spins_[i][(j+1)%L_]);
                M += spins_[i][j];
            }
        }

        // Add sample to observables
        energy_ << E/double(L_*L_);
        double m = M/double(L_*L_);
        magnetization_ << m;
        abs_magnetization_ << std::abs(M)/double(L_*L_);
        m2_ << m*m;
        m4_ << m*m*m*m;
    }

    void save(std::string const & filename){
        alps::hdf5::archive ar(filename, "wm");
        ar << alps::make_pvp("/simulation/results/"+energy_.representation(), energy_);
        ar << alps::make_pvp("/simulation/results/"+magnetization_.representation(), magnetization_);
        ar << alps::make_pvp("/simulation/results/"+abs_magnetization_.representation(), abs_magnetization_);
        ar << alps::make_pvp("/simulation/results/"+m2_.representation(), m2_);
        ar << alps::make_pvp("/simulation/results/"+m4_.representation(), m4_);
        ar << alps::make_pvp("/parameters/L", L_);
        ar << alps::make_pvp("/parameters/BETA", beta_);
        ar << alps::make_pvp("/parameters/SWEEPS", sweeps_);
        ar << alps::make_pvp("/parameters/THERMALIZATION", thermalization_);
    }

    protected:
    // Random int from the interval [0,max)
    int randint(int max) const
    {
        return static_cast<int>(max * rng_());
    }

private:
    typedef boost::mt19937 engine_type;
    typedef boost::uniform_real<> distribution_type;
    typedef boost::variate_generator<engine_type&, distribution_type> rng_type;
    engine_type eng_;
    distribution_type dist_;
    mutable rng_type rng_;

    size_t L_;
    double beta_;
    size_t sweeps_;
    size_t thermalization_;
    boost::multi_array<int,2> spins_;
    std::map< int, double > exp_table_;

    alps::RealObservable energy_;
    alps::RealObservable magnetization_;
    alps::RealObservable abs_magnetization_;
    alps::RealObservable m2_;
    alps::RealObservable m4_;

    std::string filename_;
};


class Hdf5Ising : public ::testing::TestWithParam<int> {};
TEST_P(Hdf5Ising, SimulationArchivePreservesObservablesAndParameters) {
    alps::testing::TemporaryDirectory temporary;
    constexpr std::size_t lattice = 16;
    constexpr std::size_t sweeps = 5000;
    const double beta = GetParam() / 10.;
    Simulation simulation(beta, lattice, (temporary.path() / "ising.h5").string());
    simulation.run(sweeps / 2, sweeps);
    simulation.verify_archive();
}
INSTANTIATE_TEST_SUITE_P(InverseTemperatures, Hdf5Ising, ::testing::Range(0, 11));
