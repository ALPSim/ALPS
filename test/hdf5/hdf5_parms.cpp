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

#include <iostream>
#include <boost/filesystem.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/parameter.h>

using namespace std;

int main ()
{
    
    if (boost::filesystem::exists("parms.h5") && boost::filesystem::is_regular_file("parms.h5"))
        boost::filesystem::remove("parms.h5");
    
    alps::Parameters p, p2;
    p["a"] = 10;
    p["b"] = "test";
    p["c"] = 10.;
    p["d"] = 5.;
    
    {
        alps::hdf5::archive ar("parms.h5", "a");
        ar << alps::make_pvp("/parameters", p);
    }
    {
        alps::hdf5::archive ar("parms.h5", "r");
        alps::Parameters pin;
        ar >> alps::make_pvp("/parameters", pin);
        cout << "Reading 1:" << endl << pin;
    }
    
    // "a" is modified from int to double
    // "c" is modified from double to double (but with decimals)
    p2["a"] = 10.5;
    p2["c"] = 5.2;
    {
        alps::hdf5::archive ar("parms.h5", "a");
        ar << alps::make_pvp("/parameters", p2);
    }
    {
        alps::hdf5::archive ar("parms.h5", "r");
        alps::Parameters pin;
        ar >> alps::make_pvp("/parameters", pin);
        cout << "Reading 2:" << endl << pin;
    }
    
    // "d" is modified from double to string
    p2["d"] = "newtype";
    {
        alps::hdf5::archive ar("parms.h5", "a");
        ar << alps::make_pvp("/parameters", p2);
    }
    {
        alps::hdf5::archive ar("parms.h5", "r");
        alps::Parameters pin;
        ar >> alps::make_pvp("/parameters", pin);
        cout << "Reading 3:" << endl << pin;
    }
    
    return 0;
}
