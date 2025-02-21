/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2012 by Jan Gukelberger <gukel@comp-phys.org>                     *
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
