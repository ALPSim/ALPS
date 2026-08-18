/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2012 by Jan Gukelberger <gukel@comp-phys.org>                     *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_SCHEDULER_PROGRESS_CALLBACK_HPP
#define ALPS_NGS_SCHEDULER_PROGRESS_CALLBACK_HPP

#include <alps/check_schedule.hpp>
#include <iostream>

namespace alps {
    
    class ALPS_DECL progress_callback
    {
    public:
        progress_callback(double tmin=60., double tmax=900.) : schedule_(tmin,tmax) {}
    
        /// print current progress fraction to cout if schedule says we should
        void operator()(double fraction)
        {
            if( schedule_.pending() )
            {
                schedule_.update(fraction);

                std::streamsize oldprec = std::cout.precision(3);
                std::cout << "Completed " << 100*fraction << "%.";
                if( fraction < 1. ) 
                    std::cout << " Next check in " << schedule_.check_interval() << " seconds." << std::endl;
                else
                    std::cout << " Done." << std::endl;
                std::cout.precision(oldprec);
            }
        }
    
    private:
        check_schedule schedule_;
    };

}

#endif // !defined ALPS_NGS_SCHEDULER_PROGRESS_CALLBACK_HPP
