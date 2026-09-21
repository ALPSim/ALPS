/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/params.hpp>

int main() {
    std::string const filename = "odering";
    if (boost::filesystem::exists(boost::filesystem::path(filename)))
        boost::filesystem::remove(boost::filesystem::path(filename));
    {
        alps::params parms;
        parms["a"] = 6;
        parms["x"] = 2;
        parms["b"] = 3;
        parms["w"] = 1;
        
        for (alps::params::const_iterator it = parms.begin(); it != parms.end(); ++it)
            std::cout << it->first << " " << it->second << std::endl;

        alps::hdf5::archive oar(filename, "w");
        oar["/parameters"] << parms;
    }
    std::cout << "= = = = =" << std::endl;
    {
        alps::params parms;
        alps::hdf5::archive iar(filename, "r");
        iar["/parameters"] >> parms;

        alps::params::const_iterator it = parms.begin();
        assert((it++)->first == "a");
        assert((it++)->first == "x");
        assert((it++)->first == "b");
        assert((it++)->first == "w");
    }
    boost::filesystem::remove(boost::filesystem::path(filename));
    return 0;
}
