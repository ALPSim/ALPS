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

#ifndef ALPS_NGS_MCRESULTS_HPP
#define ALPS_NGS_MCRESULTS_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/config.hpp>
#include <alps/ngs/mcresult.hpp>

#include <map>
#include <string>

namespace alps {

    class ALPS_DECL mcresults : public std::map<std::string, mcresult> {

        public: 

            mcresult & operator[](std::string const & name);

            mcresult const & operator[](std::string const & name) const;

            bool has(std::string const & name) const;

            void insert(std::string const & name, mcresult res);

            void erase(std::string const & name);

            void save(hdf5::archive & ar) const;

            void load(hdf5::archive & ar);

            void output(std::ostream & os) const;

    };

    ALPS_DECL std::ostream & operator<<(std::ostream & os, mcresults const & results);
}

#endif
