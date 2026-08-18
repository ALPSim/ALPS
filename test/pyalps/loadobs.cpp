/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Lukas Gamper <gamperl -at- gmail.com>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/hdf5/archive.hpp>
#include <alps/utility/encode.hpp>
#include <alps/alea.h>

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>


int main() {
    alps::hdf5::archive iar("loadobs.h5");

    std::vector<std::string> list = iar.list_children("/simulation/results");
    for (std::vector<std::string>::const_iterator it = list.begin(); it != list.end(); ++it) {
        iar.set_context("/simulation/results/" + iar.encode_segment(*it));
        if (iar.is_scalar("/simulation/results/" + iar.encode_segment(*it) + "/mean/value")) {
            alps::alea::mcdata<double> obs;
            obs.load(iar);
            std::cout << *it << " " << obs << std::endl;
        } else {
            alps::alea::mcdata<std::vector<double> > obs;
            obs.load(iar);
            std::cout << *it << " " << obs << std::endl;
        }
    }
}
