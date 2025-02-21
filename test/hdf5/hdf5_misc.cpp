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

#include <alps/hdf5.hpp>

#include <boost/filesystem.hpp>

class my_class {
    public:
        my_class(double v = 0): d(v) {}
        void save(alps::hdf5::archive & ar) const {
            using alps::make_pvp;
            ar << make_pvp("value", d);
        }
        void load(alps::hdf5::archive & ar) { 
            using alps::make_pvp;
            ar >> make_pvp("value", d); 
        }
    private:
        double d;
};

int main () {

    if (boost::filesystem::exists(boost::filesystem::path("data.h5")))
        boost::filesystem::remove(boost::filesystem::path("data.h5"));

    {
        alps::hdf5::archive ar("data.h5", "w");
        ar << alps::make_pvp("/value", 42);
    }
    
    {
        alps::hdf5::archive ar("data.h5");
        int i;
        ar >> alps::make_pvp("/value", i);
    }

    {
        alps::hdf5::archive ar("data.h5");
        std::string s;
        ar >> alps::make_pvp("/value", s);
    }

    {
        alps::hdf5::archive ar("data.h5", "w");
        std::vector<double> vec(5, 42);
        ar << alps::make_pvp("/path/2/vec", vec);
    }
    
    {
        std::vector<double> vec;
        // fill the vector
        alps::hdf5::archive ar("data.h5");
        ar >> alps::make_pvp("/path/2/vec", vec);
    }

    {
        std::string str("foobar");
        alps::hdf5::archive ar("data.h5", "w");
        ar << alps::make_pvp("/foo/bar", str);
    }
    
    {
        alps::hdf5::archive ar("data.h5");
        std::string str;
        ar >> alps::make_pvp("/foo/bar", str);
    }

    {
        long *d = new long[17];
        // fill the array
        alps::hdf5::archive ar("data.h5", "w");
        ar << alps::make_pvp("/c/array", d, 17);
        delete[] d;
    }

    {
        alps::hdf5::archive ar("data.h5");
        std::size_t size = ar.extent("/c/array")[0];
        long *d = new long[size];
        ar >> alps::make_pvp("/c/array", d, size);
        delete[] d;
    }

    {
        {
                my_class c(42);
                alps::hdf5::archive ar("data.h5", "w");
                ar << alps::make_pvp("/my/class", c);
        }
        {
                my_class c;
                alps::hdf5::archive ar("data.h5");
                ar >> alps::make_pvp("/my/class", c);
        }
    }

    {
        alps::hdf5::archive ar("data.h5", "w"); 
        // the parent of an attribute must exist
        ar.create_group("/foo");
        ar << alps::make_pvp("/foo/@bar", std::string("hello"));
    }

    {
        alps::hdf5::archive ar("data.h5");
        std::string str;
        ar >> alps::make_pvp("/foo/@bar", str);
    }

    boost::filesystem::remove(boost::filesystem::path("data.h5"));
    return 0;
}
