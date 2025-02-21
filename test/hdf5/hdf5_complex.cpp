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

#include <complex>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <vector>
#include <iostream>
#include <algorithm>

template<class T>
std::ostream& operator<<(std::ostream& os, std::vector<T> const & v)
{
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
    return os;
}


struct foo {
  
    std::complex<double> scalar;
    std::vector<std::complex<double> > vec;
    
    void load(alps::hdf5::archive & ar)
    {
        ar >> alps::make_pvp("scalar", scalar);
        ar >> alps::make_pvp("vector", vec);
    }
    void save(alps::hdf5::archive & ar) const
    {
        ar << alps::make_pvp("scalar", scalar);
        ar << alps::make_pvp("vector", vec);
    }
    
};
int main () {
    
    foo b;
    b.scalar = std::complex<double>(3,4);
    b.vec = std::vector<std::complex<double> >(5, std::complex<double>(0,7));
    {
        alps::hdf5::archive ar("test_hdf5_complex.h5", "w");
        ar << alps::make_pvp("/test/foo", b);
    }
    
    // check
    {
        foo t_b;
        alps::hdf5::archive ar("test_hdf5_complex.h5", "r");
        ar >> alps::make_pvp("/test/foo", t_b);
        std::cout << "scalar (write): " << b.scalar << std::endl;
        std::cout << "scalar (read): " << t_b.scalar << std::endl;
        std::cout << "vector (write): " << b.vec << std::endl;
        std::cout << "vector (read): " << t_b.vec << std::endl;
    }
    
    return 0;
}
