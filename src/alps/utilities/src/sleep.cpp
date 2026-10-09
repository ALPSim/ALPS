/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/sleep.hpp>

#ifndef ALPS_NGS_SINGLE_THREAD

#include <boost/thread.hpp>
#include <boost/thread/xtime.hpp>

namespace alps {

    void sleep(std::size_t nanoseconds) {
        // TODO: check if boost::this_thread::sleep is nicer than xtime
        boost::xtime xt;
#if BOOST_VERSION < 105000
        boost::xtime_get(&xt, boost::TIME_UTC);
#else
        boost::xtime_get(&xt, boost::TIME_UTC_);
#endif
        xt.nsec += nanoseconds;
        boost::thread::sleep(xt);
    }
}

#else

#include <ctime>
#include <stdexcept>

namespace alps {

    void sleep(std::size_t nanoseconds) {

        struct timespec tim, tim2;
        tim.tv_nsec = nanoseconds;

        if(nanosleep(&tim , &tim2) < 0)
            throw std::runtime_error("Nano sleep failed");
    }
}

#endif
