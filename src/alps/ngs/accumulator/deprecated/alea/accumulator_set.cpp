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
#include <alps/ngs/alea/accumulator_set.hpp>

XXX

namespace alps {
    namespace accumulator {

        detail::accumulator_wrapper & accumulator_set::operator[](std::string const & name) {
            if (!has(name))
                throw std::out_of_range("No observable found with the name: " + name + ALPS_STACKTRACE);
            return *(storage.find(name)->second);
        }

        detail::accumulator_wrapper const & accumulator_set::operator[](std::string const & name) const {
            if (!has(name))
                throw std::out_of_range("No observable found with the name: " + name + ALPS_STACKTRACE);
            return *(storage.find(name)->second);
        }

        bool accumulator_set::has(std::string const & name) const{
            return storage.find(name) != storage.end();
        }
        
        void accumulator_set::insert(std::string const & name, boost::shared_ptr<alps::accumulator::detail::accumulator_wrapper> ptr){
            if (has(name))
                throw std::out_of_range("There exists alrady an accumulator with the name: " + name + ALPS_STACKTRACE);
            storage.insert(make_pair(name, ptr));
        }

        void accumulator_set::save(hdf5::archive & ar) const {
            for(const_iterator it = begin(); it != end(); ++it)
                ar[it->first] = *(it->second);
        }

        void accumulator_set::load(hdf5::archive & ar) {}

        void accumulator_set::merge(accumulator_set const &) {}

        void accumulator_set::print(std::ostream & os) const {}

        void accumulator_set::reset(bool equilibrated) {
            for(iterator it = begin(); it != end(); ++it)
                it->second->reset();
        }
        
        //~ map operations
        accumulator_set::iterator accumulator_set::begin() {
            return storage.begin();
        }

        accumulator_set::iterator accumulator_set::end() {
            return storage.end();
        }

        accumulator_set::const_iterator accumulator_set::begin() const {
            return storage.begin();
        }
        
        accumulator_set::const_iterator accumulator_set::end() const {
            return storage.end();
        }
        
        void accumulator_set::clear() {
            storage.clear(); //should be ok b/c shared_ptr
        }
    }
}
