/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Mario Koenz <mkoenz@ethz.ch>                       *
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
#include <alps/ngs/alea/result_set.hpp>

XXX

namespace alps {
    namespace accumulator {

        detail::result_wrapper & result_set::operator[](std::string const & name) {
            if (!has(name))
                throw std::out_of_range("No observable found with the name: " + name + ALPS_STACKTRACE);
            return *(storage.find(name)->second);
        }

        detail::result_wrapper const & result_set::operator[](std::string const & name) const {
            if (!has(name))
                throw std::out_of_range("No observable found with the name: " + name + ALPS_STACKTRACE);
            return *(storage.find(name)->second);
        }

        bool result_set::has(std::string const & name) const {
            return storage.find(name) != storage.end();
        }

        void result_set::insert(std::string const & name, boost::shared_ptr<detail::result_wrapper> ptr) {
            if (has(name))
                throw std::out_of_range("There exists alrady a result with the name: " + name + ALPS_STACKTRACE);
            storage.insert(make_pair(name, ptr));
        }

        void result_set::save(hdf5::archive & ar) const {}

        void result_set::load(hdf5::archive & ar) {}

        void result_set::merge(result_set const &) {}

        void result_set::print(std::ostream & os) const {
            for (const_iterator it = begin(); it != end(); ++it)
                os << it->first << ": " << *(it->second);
        }

        // map operations
        result_set::iterator result_set::begin() {
            return storage.begin();
        }

        result_set::iterator result_set::end() {
            return storage.end();
        }

        result_set::const_iterator result_set::begin() const {
            return storage.begin();
        }

        result_set::const_iterator result_set::end() const {
            return storage.end();
        }

        void result_set::clear() {
            storage.clear(); //should be ok b/c shared_ptr
        }
    }
}
