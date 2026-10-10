// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
// Runtime regressions adapted from ALPS PR 166 at 519d4bd0d.
// Keep these checks independent of the native test-framework migration.
#include <alps/alea.h>
#include <alps/alea/mcdata.hpp>
#include <alps/check_schedule.hpp>
#include <alps/parapack/integer_range.h>
#include <limits>
#include <alps/numeric/vector_valarray_conversion.hpp>
#include <alps/utility/temporary_filename.hpp>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#ifndef _WIN32
#include <sys/stat.h>
#endif

void require(bool valid, const char* message) {
    if (!valid) throw std::runtime_error(message);
}

int main() {
    const alps::integer_range<int> full((std::numeric_limits<int>::min)(),
                                        (std::numeric_limits<int>::max)());
    require(!full.empty() && unify(full, full) == full, "full-width range union");
    bool overflow = false;
    try { (void)full.size(); } catch (std::overflow_error const&) { overflow = true; }
    require(overflow, "unrepresentable range size");
    alps::check_schedule schedule;
    require(schedule.check_interval() == 0., "initial check interval");
    require(schedule.pending(), "first check must be pending");
    schedule.update(0.);
    require(schedule.check_interval() == 60., "first update interval");
    require(alps::numeric::vector2valarray(std::vector<double>{}).size() == 0,
            "empty vector conversion");
    require(alps::numeric::valarray2vector(std::valarray<double>{}).empty(),
            "empty valarray conversion");

    alps::RealObservable observable("scalar");
    for (int i = 0; i < 1024; ++i) observable << (1. + i % 10);
    const alps::alea::mcdata<double> original(observable);
    const double mean = original.mean();
    const auto negative = -original;
    require(negative.mean() == -mean, "negation mean");
    require(original.mean() == mean, "negation must not mutate source");
    require(negative.error() >= 0., "negation uncertainty");
    auto reciprocal = original;
    reciprocal.divide(-2.);
    require(reciprocal.error() >= 0., "reciprocal uncertainty");
    alps::RealObsevaluator evaluator(observable);
    evaluator *= -2.;
    require(evaluator.error() >= 0., "negative scaling uncertainty");
    evaluator /= -2.;
    require(evaluator.error() >= 0., "negative division uncertainty");

    alps::ObservableSet first, second;
    first << alps::IntHistogramObservable("histogram", 0, 10);
    second << alps::IntHistogramObservable("histogram", 0, 10);
    first << second;
    const alps::IntHistogramObsevaluator histogram(first["histogram"]);
    require(histogram.count() == 0 && histogram.size() == 10, "empty histogram range");
    for (std::size_t i = 0; i < histogram.size(); ++i)
        require(histogram[i] == 0, "empty histogram bin");
    first.clear();
    require(first.empty(), "observable clear");
#ifndef _WIN32
    const mode_t previous = umask(0);
    const auto filename = alps::temporary_filename("alps-runtime-regression-");
    umask(previous);
    struct stat information{};
    const int status = stat(filename.c_str(), &information);
    std::remove(filename.c_str());
    require(status == 0 && (information.st_mode & 0777) == 0600, "private temporary file");
#endif
}
