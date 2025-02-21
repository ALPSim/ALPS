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
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>

#include <boost/filesystem.hpp>

#include <iostream>
using namespace std;

int main()
{
    if (boost::filesystem::exists(boost::filesystem::path("test_hdf5_memory.h5")))
        boost::filesystem::remove(boost::filesystem::path("test_hdf5_memory.h5"));
    {
        alps::hdf5::archive oa("test_hdf5_memory.h5", "w");
        std::vector<std::complex<double> > foo(3);
        std::vector<double> foo2(3);
        oa << alps::make_pvp("/foo", foo);
        oa << alps::make_pvp("/foo2", foo2);
    }
    
    {
 
        std::vector<double> foo, foo2;
        try {
			alps::hdf5::archive ia("test_hdf5_memory.h5");
            ia >> alps::make_pvp("/foo", foo);
            ia >> alps::make_pvp("/foo2", foo2);
        } catch (exception e) {
            cout << "Exception caught: no complex value" << endl;
            boost::filesystem::remove(boost::filesystem::path("test_hdf5_memory.h5"));
            return EXIT_SUCCESS;
        }
    }
    boost::filesystem::remove(boost::filesystem::path("test_hdf5_memory.h5"));
    return EXIT_SUCCESS;
}