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
 
#include <boost/filesystem.hpp>

#include <string>
#include <iostream>

typedef enum { PLUS, MINUS } enum_type;

void save(
      alps::hdf5::archive & ar
    , std::string const & path
    , enum_type const & value
    , std::vector<std::size_t> size = std::vector<std::size_t>()
    , std::vector<std::size_t> chunk = std::vector<std::size_t>()
    , std::vector<std::size_t> offset = std::vector<std::size_t>()
) {
    switch (value) {
        case PLUS: ar << alps::make_pvp(path, std::string("plus")); break;
        case MINUS: ar << alps::make_pvp(path, std::string("minus")); break;
    }
}

void load(
      alps::hdf5::archive & ar
    , std::string const & path
    , enum_type & value
    , std::vector<std::size_t> chunk = std::vector<std::size_t>()
    , std::vector<std::size_t> offset = std::vector<std::size_t>()
) {
    std::string s;
    ar >> alps::make_pvp(path, s);
    value = (s == "plus" ? PLUS : MINUS);
}

int main() {

    std::string const filename = "example_enum_as_class.h5";

    if (boost::filesystem::exists(boost::filesystem::path(filename)))
        boost::filesystem::remove(boost::filesystem::path(filename));

    enum_type read, write = PLUS;

    {
        alps::hdf5::archive ar(filename, "a");
        ar << alps::make_pvp("/enum", write);
    }

    {
        alps::hdf5::archive ar(filename);
        ar >> alps::make_pvp("/enum", read);
    }

    boost::filesystem::remove(boost::filesystem::path(filename));

    bool match = (read == write);
    std::cout << (match ? "SUCCESS" : "FAILURE") << std::endl;
    return match ? EXIT_SUCCESS : EXIT_FAILURE;
}
