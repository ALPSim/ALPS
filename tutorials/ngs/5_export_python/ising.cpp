// Copyright (C) 2010-2012 by Lukas Gamper
//               2026      by the ALPS collaboration
// SPDX-License-Identifier: MIT

#include "ising.hpp"

#include <alps/hdf5/archive.hpp>
#include <alps/ngs.hpp>

#include <algorithm>

ising_sim::ising_sim(parameters_type const & parameters,
                     std::size_t seed_offset)
    : alps::mcbase(parameters, seed_offset),
      total_sweeps_(parameters["SWEEPS"] | 10) {
    measurements << alps::accumulator::RealObservable("Magnetization");
}

void ising_sim::update() {
    state_ = random() < 0.5 ? -1.0 : 1.0;
    ++sweeps_;
}

void ising_sim::measure() {
    measurements["Magnetization"] << state_;
}

double ising_sim::fraction_completed() const {
    return std::min(1.0, static_cast<double>(sweeps_) / total_sweeps_);
}

void ising_sim::save(alps::hdf5::archive & archive) const {
    alps::mcbase::save(archive);
    archive["/checkpoint/sweeps"] << sweeps_;
    archive["/checkpoint/state"] << state_;
}

void ising_sim::load(alps::hdf5::archive & archive) {
    alps::mcbase::load(archive);
    archive["/checkpoint/sweeps"] >> sweeps_;
    archive["/checkpoint/state"] >> state_;
}
