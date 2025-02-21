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
