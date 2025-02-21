/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
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

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/complex.hpp>

#include <boost/filesystem.hpp>

#include <iostream>

using namespace alps;

int main() {

    std::string const filename = "ngs.h5";
    if (boost::filesystem::exists(boost::filesystem::path(filename)))
        boost::filesystem::remove(boost::filesystem::path(filename));
    {
        hdf5::archive ar(filename, "a");
        ar << make_pvp("/to/to", 3.14159);
    }
    {
        hdf5::archive ar(filename, "r");
        double value;
        ar >> make_pvp("/to/to", value);
        std::cout << value << std::endl;
    }

    {
        hdf5::archive ar(filename, "a");
        ar << make_pvp("/to/my/vec/in/a/very/deep/path", std::vector<double>(17, 15.141));
    }
    {
        hdf5::archive ar(filename, "r");
        std::vector<unsigned> value;
        ar >> make_pvp("/to/my/vec/in/a/very/deep/path", value);
        std::cout << value[0] << std::endl;
    }

    {
        hdf5::archive ar(filename, "a");
        ar << make_pvp("/to/to", std::complex<double>(3.14159, 12.34));
    }
    {
        hdf5::archive ar(filename, "r");
        std::complex<double> value;
        ar >> make_pvp("/to/to", value);
        std::cout << value.real() << " " << value.imag() << std::endl;
    }

    {
        hdf5::archive ar(filename, "a");
        ar << make_pvp("/to/str", std::string("asdf"));
    }
    {
        hdf5::archive ar(filename, "r");
        std::string value;
        ar >> make_pvp("/to/str", value);
        std::cout << value << std::endl;
    }

    {
        hdf5::archive ar(filename, "a");
        ar << make_pvp("/to/char", "asdf");
    }
    {
        hdf5::archive ar(filename, "r");
        std::string value;
        ar >> make_pvp("/to/char", value);
        std::cout << value << std::endl;
    }
    {
        hdf5::archive ar(filename, "r");
        std::cout << (ar.is_datatype<double>("/to/to") ? "true" : "false") << std::endl;
    }
    boost::filesystem::remove(boost::filesystem::path(filename));
    return 0;
}
