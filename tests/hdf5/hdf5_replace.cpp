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
#include <alps/hdf5/vector.hpp>

#include <vector>

using namespace std;

int main () {

   vector<double> vec(100, 10.);

   {
       alps::hdf5::archive ar("res.h5", "a");
       ar << alps::make_pvp("/vec2", vec);
   }
   {
       alps::hdf5::archive ar("res.h5", "w");
       ar << alps::make_pvp("/vec", vec);
   }
   {
       vector<double> tmp;
       alps::hdf5::archive ar("res.h5");
       ar >> alps::make_pvp("/vec2", tmp);
   }

    return 0;
}
