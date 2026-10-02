/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>,                  *
 *                              Synge Todo <wistaria@comp-phys.org>                *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
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
