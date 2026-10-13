/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once
#include <alps/ngs.hpp>
#include <alps/mcbase.hpp>
#include <cmath>
#include <vector>

// Retain the original exp(-x*x) scalar/vector simulation, with defined state.
// Record the sample sum independently of the observable implementation.
class sum_simulation : public alps::mcbase {
public:
    sum_simulation(parameters_type const& parameters, std::size_t seed_offset = 42)
        : alps::mcbase(parameters, seed_offset), total_(parameters["COUNT"]) {
        measurements << alps::accumulator::RealObservable("SValue")
                     << alps::accumulator::RealVectorObservable("VValue");
    }
    void update() override {
        const double x = random();
        value_ = std::exp(-x * x);
    }
    void measure() override {
        ++count;
        sum += value_;
        measurements["SValue"] << value_;
        measurements["VValue"] << std::vector<double>(3, value_);
    }
    double fraction_completed() const override { return double(count) / total_; }
    int count = 0;
    double sum = 0;
private:
    int total_;
    double value_ = 0;
};
