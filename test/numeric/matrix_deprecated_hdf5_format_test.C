/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013        by Michele Dolfi <dolfim@phys.ethz.ch>,               *
 *                              Andreas Hehn <hehn@phys.ethz.ch>                   *
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
#include <alps/numeric/matrix.hpp>
#include <alps/version.h>
#include <alps/hdf5.hpp>


int main() {
    boost::filesystem::path   infile(ALPS_SRCDIR);
    infile = infile / "test" / "numeric" / "matrix_deprecated_hdf5_format_test.h5";
    if (!boost::filesystem::exists(infile))
    {
        std::cout << "Reference file " << infile << " not found." << std::endl;
        return -1;
    }
    alps::numeric::matrix<double> m;
    alps::hdf5::archive ar(infile.native(), "r");
    ar["/matrix_old_hdf5_format"] >> m;
    std::cout << "Matrix " << num_rows(m) << "x" << num_cols(m) << ":\n";
    std::cout << "capacity: " <<m.capacity().first << "x" << m.capacity().second <<"\n";
    std::cout << "data:\n" << m;
    return 0;
}
