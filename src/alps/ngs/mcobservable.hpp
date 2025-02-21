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

#ifndef ALPS_NGS_MCOBSERVABLE_HPP
#define ALPS_NGS_MCOBSERVABLE_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/config.hpp>

#include <alps/alea/observable_fwd.hpp>

#include <map>
#include <iostream>

namespace alps {

    class ALPS_DECL mcobservable {

        public:

            mcobservable();
            mcobservable(Observable const * obs);
            mcobservable(mcobservable const & rhs);

            virtual ~mcobservable();

            mcobservable & operator=(mcobservable rhs);

            Observable * get_impl();

            Observable const * get_impl() const;

            std::string const & name() const;

            template<typename T> mcobservable & operator<<(T const & value);

            void save(hdf5::archive & ar) const;
            void load(hdf5::archive & ar);

            void merge(mcobservable const &);

            void output(std::ostream & os) const;

        private:

            Observable * impl_;
            static std::map<Observable *, std::size_t> ref_cnt_;

    };

    std::ostream & operator<<(std::ostream & os, mcobservable const & obs);

}

#endif
