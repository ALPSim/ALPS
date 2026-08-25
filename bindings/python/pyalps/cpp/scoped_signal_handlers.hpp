// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT

#ifndef PYALPS_SCOPED_SIGNAL_HANDLERS_HPP
#define PYALPS_SCOPED_SIGNAL_HANDLERS_HPP

#include <alps/ngs/signal.hpp>

#if !defined(BOOST_MSVC) && !defined(ALPS_NGS_NO_SIGNALS)
#include <array>
#include <csignal>
#include <cstring>
#endif

namespace pyalps {

// ALPS applications own process signal handlers while they run. Python is an
// embedded runtime, however, so its handlers must be put back before control
// returns to the interpreter. Reinstalling the ALPS handlers here also makes
// repeated solve() calls work after the preceding guard restored Python's.
class scoped_signal_handlers {
public:
    scoped_signal_handlers() {
#if !defined(BOOST_MSVC) && !defined(ALPS_NGS_NO_SIGNALS)
        for (std::size_t i = 0; i < signal_numbers_.size(); ++i)
            saved_[i].valid = sigaction(signal_numbers_[i], NULL, &saved_[i].action) == 0;

        struct sigaction action;
        std::memset(&action, 0, sizeof(action));
        action.sa_handler = &alps::ngs::signal::slot;
        for (std::size_t i = 0; i < termination_signal_count_; ++i)
            sigaction(signal_numbers_[i], &action, NULL);

        action.sa_handler = &alps::ngs::signal::segfault;
        for (std::size_t i = termination_signal_count_; i < signal_numbers_.size(); ++i)
            sigaction(signal_numbers_[i], &action, NULL);
#endif
    }

    ~scoped_signal_handlers() {
#if !defined(BOOST_MSVC) && !defined(ALPS_NGS_NO_SIGNALS)
        for (std::size_t i = 0; i < signal_numbers_.size(); ++i)
            if (saved_[i].valid)
                sigaction(signal_numbers_[i], &saved_[i].action, NULL);
#endif
    }

    scoped_signal_handlers(scoped_signal_handlers const &) = delete;
    scoped_signal_handlers & operator=(scoped_signal_handlers const &) = delete;

private:
#if !defined(BOOST_MSVC) && !defined(ALPS_NGS_NO_SIGNALS)
    struct saved_action {
        struct sigaction action;
        bool valid = false;
    };

    static constexpr std::array<int, 8> signal_numbers_ = {{
        SIGINT, SIGTERM, SIGXCPU, SIGQUIT, SIGUSR1, SIGUSR2, SIGSEGV, SIGBUS
    }};
    static constexpr std::size_t termination_signal_count_ = 6;
    std::array<saved_action, signal_numbers_.size()> saved_;
#endif
};

#if !defined(BOOST_MSVC) && !defined(ALPS_NGS_NO_SIGNALS)
constexpr std::array<int, 8> scoped_signal_handlers::signal_numbers_;
constexpr std::size_t scoped_signal_handlers::termination_signal_count_;
#endif

} // namespace pyalps

#endif
