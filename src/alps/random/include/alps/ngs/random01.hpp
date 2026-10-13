/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
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
