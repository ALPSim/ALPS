/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>,                  *
 *                              Synge Todo <wistaria@comp-phys.org>                *
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
#include <alps/ngs/boost_mpi.hpp>
#include <alps/stop_callback.hpp>

namespace alps {

    stop_callback::stop_callback(std::size_t timelimit)
        : limit(timelimit)
        , start(boost::chrono::high_resolution_clock::now())
    {}

#ifdef ALPS_HAVE_MPI
    stop_callback::stop_callback(boost::mpi::communicator const & cm, std::size_t timelimit)
        : limit(timelimit), start(boost::chrono::high_resolution_clock::now()), comm(cm)
    {}
#endif

    bool stop_callback::operator()() {
#ifdef ALPS_HAVE_MPI
        if (comm) {
            bool to_stop;
            if (comm->rank() == 0)
                to_stop = !signals.empty() || (limit.count() > 0 && boost::chrono::high_resolution_clock::now() > start + limit);
            broadcast(*comm, to_stop, 0);
            return to_stop;
        } else
#endif
            return !signals.empty() || (limit.count() > 0 && boost::chrono::high_resolution_clock::now() > start + limit);
    }

#ifdef ALPS_HAVE_MPI
    stop_callback_mpi::stop_callback_mpi(boost::mpi::communicator const & cm, std::size_t timelimit)
        : comm(cm), limit(timelimit), start(boost::chrono::high_resolution_clock::now())
    {}

    bool stop_callback_mpi::operator()() {
        bool to_stop;
        if (comm.rank() == 0)
            to_stop = !signals.empty() 
               || (limit.count() > 0 && boost::chrono::high_resolution_clock::now() > start + limit);
        broadcast(comm, to_stop, 0);
        return to_stop;
    }
#endif
}
