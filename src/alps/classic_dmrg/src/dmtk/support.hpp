// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Standard-library replacements for the ALPS utility helpers that the alpsx
// dmtk used (alps::Stopwatch, alps::format_elapsed and
// alps::temporary_filename), so the component links only ALPS::models.
#pragma once

#include <atomic>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

namespace dmtk {

// Mirrors the alpsx Stopwatch: restart() resets both the start and the mark,
// mark() only the mark; elapsed() and lap_elapsed() report seconds since them.
class Stopwatch {
public:
    Stopwatch() : start_(clock::now()), mark_(start_) {}
    void restart() { start_ = mark_ = clock::now(); }
    void mark() { mark_ = clock::now(); }
    double elapsed() const { return seconds(clock::now() - start_); }
    double lap_elapsed() const { return seconds(clock::now() - mark_); }

private:
    using clock = std::chrono::steady_clock;
    static double seconds(clock::duration d) { return std::chrono::duration<double>(d).count(); }
    clock::time_point start_, mark_;
};

inline std::string format_elapsed(double seconds) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(3) << seconds << " s";
    return out.str();
}

// Scratch files live in a directory created privately for one invocation, so
// a process-wide counter makes the names unique; nothing is created here.
inline std::string temporary_filename(std::string const& prefix) {
    static std::atomic<unsigned long> counter{0};
    return prefix + "." + std::to_string(counter++) + ".tmp";
}

} // namespace dmtk
