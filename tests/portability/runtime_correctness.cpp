// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
// Runtime regressions adapted from ALPS PR 166 at 519d4bd0d.
// Preserve the SDK follow-up regression coverage during test consolidation.
#include <gtest/gtest.h>
#include <alps/mc/check_schedule.hpp>
#include <alps/parapack/integer_range.h>
#include <limits>
#include <alps/numeric/vector_valarray_conversion.hpp>
#include <alps/utilities/temporary_filename.hpp>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#ifndef _WIN32
#include <sys/stat.h>
#endif

void require(bool valid, const char* message) {
    if (!valid) throw std::runtime_error(message);
}

TEST(RuntimeCorrectness, PreservesSDKRegressions) {
    const alps::integer_range<int> full((std::numeric_limits<int>::min)(),
                                        (std::numeric_limits<int>::max)());
    require(!full.empty() && unify(full, full) == full, "full-width range union");
    bool overflow = false;
    try { (void)full.size(); } catch (std::overflow_error const&) { overflow = true; }
    require(overflow, "unrepresentable range size");
    alps::check_schedule schedule(60., 900.);
    require(schedule.pending(), "first check must be pending");
    schedule.update(0.);
    require(!schedule.pending(), "first update schedules the next check");
    require(alps::numeric::vector2valarray(std::vector<double>{}).size() == 0,
            "empty vector conversion");
    require(alps::numeric::valarray2vector(std::valarray<double>{}).empty(),
            "empty valarray conversion");

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
