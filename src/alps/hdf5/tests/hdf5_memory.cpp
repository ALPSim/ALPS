/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
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