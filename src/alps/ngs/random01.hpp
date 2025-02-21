/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
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

#ifndef ALPS_NGS_RANDOM01_HPP
#define ALPS_NGS_RANDOM01_HPP

#include <alps/hdf5/archive.hpp>

#include <boost/random.hpp>

#include <string>
#include <sstream>

namespace alps {

    struct random01 : public boost::variate_generator<boost::mt19937, boost::uniform_01<double> > {
        random01(int seed = 42)
            : boost::variate_generator<boost::mt19937, boost::uniform_01<double> >(boost::mt19937(seed), boost::uniform_01<double>())
        {}

        void save(alps::hdf5::archive & ar) const { // TODO: move this to hdf5 archive!
            std::ostringstream os;
            os << this->engine();
            ar["engine"] << os.str();
        }

        void load(alps::hdf5::archive & ar) { // TODO: move this to hdf5 archive!
            std::string state;
            ar["engine"] >> state;
            std::istringstream is(state);
            is >> this->engine();
        }
    };

}

#endif 
