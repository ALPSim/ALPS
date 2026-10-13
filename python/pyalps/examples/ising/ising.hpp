// Copyright (C) 2010-2012 by Lukas Gamper
//               2026      by the ALPS collaboration
// SPDX-License-Identifier: MIT
#ifndef ALPS_TUTORIAL_EXPORTED_ISING_HPP
#define ALPS_TUTORIAL_EXPORTED_ISING_HPP

#include <alps/mcbase.hpp>

#include <cstddef>

class ising_sim : public alps::mcbase {
public:
    explicit ising_sim(parameters_type const & parameters,
                       std::size_t seed_offset = 0);

    void update() override;
    void measure() override;
    double fraction_completed() const override;
    void save(alps::hdf5::archive & archive) const override;
    void load(alps::hdf5::archive & archive) override;

private:
    std::size_t sweeps_ = 0;
    std::size_t total_sweeps_ = 1;
    double state_ = 1.0;
};

#endif
