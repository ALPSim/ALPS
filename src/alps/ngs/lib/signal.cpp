/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/signal.hpp>
#include <alps/ngs/stacktrace.hpp>
#include <alps/hdf5/archive.hpp>

#include <cstring>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <cstdio>

#include <signal.h>

namespace alps {
    namespace ngs {

        signal::signal() {
            #if not ( defined BOOST_MSVC || defined ALPS_NGS_NO_SIGNALS )
                static bool initialized;
                if (!initialized) {
                    initialized = true;

                    static struct sigaction action;
                    memset(&action, 0, sizeof(action));
                    action.sa_handler = &signal::slot;
                    sigaction(SIGINT, &action, NULL);
                    sigaction(SIGTERM, &action, NULL);
                    sigaction(SIGXCPU, &action, NULL);
                    sigaction(SIGQUIT, &action, NULL);
                    sigaction(SIGUSR1, &action, NULL);
                    sigaction(SIGUSR2, &action, NULL);
                    sigaction(SIGSTOP, &action, NULL);
                    sigaction(SIGKILL, &action, NULL);
                }
            #endif
            listen();
        }

        bool signal::empty() {
            return end_  == begin_;
        }

        int signal::top() {
            return signals_[(end_ - 1) & 0x1F];
        }

        void signal::pop() {
            --end_ &= 0x1F;
        }

        void signal::listen() {
            #if not ( defined BOOST_MSVC || defined ALPS_NGS_NO_SIGNALS )
                static bool initialized;
                if (!initialized) {
                    initialized = true;

                    static struct sigaction action;
                    memset(&action, 0, sizeof(action));
                    action.sa_handler = &signal::segfault;
                    sigaction(SIGSEGV, &action, NULL);
                    sigaction(SIGBUS, &action, NULL);
                }
            #endif
        }

        void signal::slot(int signal) {
            fprintf(stderr, "Received signal %i\n", signal);
            signals_[end_] = signal;
            ++end_ &= 0x1F;
            if (begin_ == end_)
                ++begin_ &= 0x1F;
        }

        void signal::segfault(int signal) {
            hdf5::archive::abort();
            fprintf(stderr, "Abort by signal %i\n", signal);
            std::cerr << ALPS_STACKTRACE;
            std::abort();
        }

        std::size_t signal::begin_ = 0;
        std::size_t signal::end_ = 0;
        boost::array<int, 32> signal::signals_;
    }
}
