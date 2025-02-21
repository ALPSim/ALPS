/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
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

#include <vector>

using namespace std;
using namespace alps;

int main () {

    for (std::size_t i = 0; i < 100; ++i)
        if (boost::filesystem::exists(boost::filesystem::path("large" + alps::cast<std::string>(i) + ".h5")))
            boost::filesystem::remove(boost::filesystem::path("large" + alps::cast<std::string>(i) + ".h5"));

    hdf5::archive ar("large%d.h5", "al");
    for (unsigned long long s = 1; s < (1ULL << 29); s <<= 1) {
        std::cout << s << std::endl;
        vector<double> vec(s, 10.);
        ar << make_pvp("/" + cast<std::string>(s), vec);
    }

    for (std::size_t i = 0; i < 100; ++i)
        if (boost::filesystem::exists(boost::filesystem::path("large" + alps::cast<std::string>(i) + ".h5")))
            boost::filesystem::remove(boost::filesystem::path("large" + alps::cast<std::string>(i) + ".h5"));
    return 0;
}
