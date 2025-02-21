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
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/multi_array.hpp>

#include <boost/multi_array.hpp>
#include <boost/filesystem.hpp>

#include <vector>

using namespace std;
using boost::multi_array;

int main()
{
    if (boost::filesystem::exists(boost::filesystem::path("test_hdf5_multi_array.h5")))
        boost::filesystem::remove(boost::filesystem::path("test_hdf5_multi_array.h5"));

    multi_array<double,2> a( boost::extents[3][3] );
    multi_array<double,2> b( boost::extents[4][4] );

    // Write
    {
        alps::hdf5::archive ar("test_hdf5_multi_array.h5","a");
        vector< multi_array<double,2> > v(2,a);
        ar << alps::make_pvp("uniform",v);
        v.push_back(b);
        ar << alps::make_pvp("nonuniform",v);
    }

    // Read
    {
        alps::hdf5::archive ar("test_hdf5_multi_array.h5","r");
        vector< multi_array<double,2> > w;
        ar >> alps::make_pvp("nonuniform",w);
        cout << "read nonuniform" << endl;
        ar >> alps::make_pvp("uniform",w); // throws runtime_error
        cout << "read uniform" << endl;
    }
    
    if (boost::filesystem::exists(boost::filesystem::path("test_hdf5_multi_array.h5")))
        boost::filesystem::remove(boost::filesystem::path("test_hdf5_multi_array.h5"));
    return 0;
}
