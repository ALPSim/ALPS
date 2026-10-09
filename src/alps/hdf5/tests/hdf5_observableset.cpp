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


#include <alps/alea.h>

int main(int argc, char** argv)
{
    try {
        
        alps::ObservableSet measurements_;
        measurements_ << alps::RealObservable("E");
        
        alps::hdf5::archive ar("test_observableset.h5", "a");
        measurements_.get<alps::RealObservable>("E") << 1;
        ar << alps::make_pvp("/simulation/results/", measurements_);
        
    } catch (std::exception & e) {
        std::cerr << e.what() << std::endl;
        exit(1);
    }
    
    return 0;
}
