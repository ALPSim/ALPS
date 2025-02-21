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

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/mcresults.hpp>
#include <alps/ngs/stacktrace.hpp>

#include <alps/alea/observableset.h>

#include <stdexcept>

namespace alps {

    mcresult & mcresults::operator[](std::string const & name) {
        if (!has(name))
            throw std::out_of_range("No result found with the name: " + name + "\n" + ALPS_STACKTRACE);
        return std::map<std::string, mcresult>::find(name)->second;
    }

    mcresult const & mcresults::operator[](std::string const & name) const {
        if (!has(name))
            throw std::out_of_range("No result found with the name: " + name + "\n" + ALPS_STACKTRACE);
        return std::map<std::string, mcresult>::find(name)->second;
    }

    bool mcresults::has(std::string const & name) const {
        return std::map<std::string, mcresult>::find(name) != std::map<std::string, mcresult>::end();
    }

    void mcresults::insert(std::string const & name, mcresult res) {
        if (has(name))
            throw std::out_of_range("There exists alrady a result with the name: " + name + "\n" + ALPS_STACKTRACE);
        std::map<std::string, mcresult>::insert(make_pair(name, res));
    }

    void mcresults::erase(std::string const & name) {
        if (!has(name))
            throw std::out_of_range("There is no result with the name: " + name + "\n" + ALPS_STACKTRACE);
        std::map<std::string, mcresult>::erase(name);
    }

    void mcresults::save(hdf5::archive & ar) const {
        for(std::map<std::string, mcresult>::const_iterator it = std::map<std::string, mcresult>::begin(); it != std::map<std::string, mcresult>::end(); ++it)
            if (it->second.count())
                ar
                    << make_pvp(ar.encode_segment(it->first), it->second)
                ;
    }

    void mcresults::load(hdf5::archive & ar)  {
        ObservableSet set;
        // TODO: do not use hard coded path!
        ar >> make_pvp("/simulation/realizations/0/clones/0/results", set);
        for(ObservableSet::const_iterator it = set.begin(); it != set.end(); ++it)
            insert(it->first, mcresult(it->second));
    }

    void mcresults::output(std::ostream & os) const {
        for(std::map<std::string, mcresult>::const_iterator it = std::map<std::string, mcresult>::begin(); it != std::map<std::string, mcresult>::end(); ++it)
            os << std::fixed << std::setprecision(5) << it->first << ": " << it->second << std::endl;
    }

    std::ostream & operator<<(std::ostream & os, mcresults const & results) {
        results.output(os);
        return os;
    }

}
