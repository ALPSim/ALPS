/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Michele Dolfi <dolfim@phys.ethz.ch>                *
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
#include <alps/utility/vectorio.hpp>

#include <boost/filesystem.hpp>

#include <vector>
#include <complex>
#include <iostream>

template <class T>
std::ostream& operator<<(std::ostream& os, std::vector<T> const& v)
{
	os << "[" <<alps::write_vector(v, " ", 6) << "]";
	return os;
}

int main() {

    if (boost::filesystem::exists("real_complex_vec.h5") && boost::filesystem::is_regular_file("real_complex_vec.h5"))
        boost::filesystem::remove("real_complex_vec.h5");

    try {
        const int size = 6;

        std::vector<double> v(size, 3.2);

        std::cout << "v: " << v << std::endl;

        {
            alps::hdf5::archive ar("real_complex_vec.h5", "w");
            ar["/vec"] << v;
        }

        std::vector<std::complex<double> > w;
        {
            alps::hdf5::archive ar("real_complex_vec.h5", "r");
            ar["/vec"] >> w;
        }

        std::cout << "w: " << w << std::endl;
        
        boost::filesystem::remove("real_complex_vec.h5");
        
		bool passed = true;
		for (int i=0; passed && i<size; ++i)
			passed = (v[i] == w[i]);

        std::cout << "Test status checked element by element." << std::endl;
        return (passed) ? EXIT_SUCCESS : EXIT_FAILURE;

    } catch (alps::hdf5::archive_error) {
        boost::filesystem::remove("real_complex_vec.h5");
        std::cout << "Test passed because Exception was thrown." << std::endl;
        return EXIT_SUCCESS;
    }

}