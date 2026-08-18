/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013        by Michele Dolfi <dolfim@phys.ethz.ch>,               *
 *                              Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
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
