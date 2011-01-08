/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                           Matthias Troyer <troyer@comp-phys.org>                *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_MCRESULTS_HPP
#define ALPS_NGS_MCRESULTS_HPP

#include <alps/ngs/mcresult.hpp>

#include <alps/config.h>
#include <alps/hdf5_fwd.hpp>

#include <map>
#include <string>

namespace alps {

    class mcresults : public std::map<std::string, mcresult> {

        public: 

            mcresult & operator[](std::string const & name);

            mcresult const & operator[](std::string const & name) const;

            bool has(std::string const & name) const;

            void insert(std::string const & name, mcresult res);

            void serialize(hdf5::iarchive & ar);

            void serialize(hdf5::oarchive & ar) const;

            void output(std::ostream & os) const;

    };

    std::ostream & operator<<(std::ostream & os, mcresults const & mcresults);
}

#endif
