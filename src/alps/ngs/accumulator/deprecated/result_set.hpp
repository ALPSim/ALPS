/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013 by Lukas Gamper <gamperl@gmail.ch>                           *
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

#ifndef ALPS_NGS_ALEA_RESULT_SET_HPP
#define ALPS_NGS_ALEA_RESULT_SET_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/alea/wrapper/result_wrapper.hpp>

#include <map>
#include <string>

namespace alps {
    namespace accumulator {

        class ALPS_DECL result_set {

            public: 
                typedef std::map<std::string, boost::shared_ptr<detail::result_wrapper> > map_type;

                typedef map_type::iterator iterator;
                typedef map_type::const_iterator const_iterator;

                detail::result_wrapper & operator[](std::string const & name);

                detail::result_wrapper const & operator[](std::string const & name) const;

                bool has(std::string const & name) const;

                void insert(std::string const & name, boost::shared_ptr<detail::result_wrapper> ptr);

                void save(hdf5::archive & ar) const;

                void load(hdf5::archive & ar);

                void merge(result_set const &);

                void print(std::ostream & os) const;

                iterator begin();
                iterator end();
                const_iterator begin() const;
                const_iterator end() const;
                void clear();

            private:
                map_type storage;
        };

        inline std::ostream & operator<<(std::ostream & out, result_set const & arg) {
            arg.print(out);
            return out;
        }

    }
}

#endif
