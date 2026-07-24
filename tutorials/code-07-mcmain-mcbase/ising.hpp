/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_TUTORIAL_ISING_HPP
#define ALPS_TUTORIAL_ISING_HPP

#include <alps/mcbase.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>

#include <boost/function.hpp>
#include <boost/filesystem/path.hpp>

#include <vector>
#include <string>

class ALPS_DECL ising_sim : public alps::mcbase {

    public:
        
        ising_sim(parameters_type const & parms, std::size_t seed_offset = 0);

        virtual void update();
        virtual void measure();
        virtual double fraction_completed() const;

        using alps::mcbase::save;
        virtual void save(alps::hdf5::archive & ar) const;

        using alps::mcbase::load;
        virtual void load(alps::hdf5::archive & ar);

    private:
        
        int length;
        int sweeps;
        int thermalization_sweeps;
        int total_sweeps;
        double beta;
        std::vector<int> spins;
};

#endif
