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

#include <chrono>
#include <iostream>

namespace alps {

    /// Prints the completed fraction to std::cout, at most as often as its
    /// schedule allows. An ALPS addition beside ALPSCore's mc module.
    ///
    /// The schedule is the one ALPS's alps::check_schedule had before
    /// ALPSCore's replaced it: ALPSCore's cannot report the interval to the
    /// next check, which this prints.
    class progress_callback
    {
    public:
        progress_callback(double tmin=60., double tmax=900.)
        :   min_check_(tmin)
        ,   max_check_(tmax)
        ,   start_fraction_(0.)
        ,   next_check_(duration::zero())
        {}

        /// print current progress fraction to cout if schedule says we should
        void operator()(double fraction)
        {
            if( pending() )
            {
                update(fraction);

                std::streamsize oldprec = std::cout.precision(3);
                std::cout << "Completed " << 100*fraction << "%.";
                if( fraction < 1. )
                    std::cout << " Next check in " << next_check_.count() << " seconds." << std::endl;
                else
                    std::cout << " Done." << std::endl;
                std::cout.precision(oldprec);
            }
        }

    private:
        typedef std::chrono::steady_clock clock;
        typedef clock::time_point time_point;
        typedef std::chrono::duration<double> duration;

        bool pending() const
        {
            return clock::now() > last_check_time_ + next_check_;
        }

        void update(double fraction)
        {
            time_point now = clock::now();

            // first check
            if( start_time_ == time_point() )
            {
                start_time_ = now;
                start_fraction_ = fraction;
            }

            if( fraction > start_fraction_ )
            {
                // estimate remaining time; propose to run 1/4 of that time
                duration old_check = next_check_;
                next_check_ = 0.25 * (1 - fraction) * (now - start_time_) / (fraction - start_fraction_);
                if( next_check_ > 2*old_check ) next_check_ = 2 * old_check;
                if( next_check_ < min_check_ ) next_check_ = min_check_;
                if( next_check_ > max_check_ ) next_check_ = max_check_;
            }
            else
                next_check_ = min_check_;

            last_check_time_ = now;
        }

        duration min_check_;
        duration max_check_;

        time_point start_time_;
        time_point last_check_time_;
        double start_fraction_;
        duration next_check_;
    };

}

#endif // !defined ALPS_NGS_SCHEDULER_PROGRESS_CALLBACK_HPP
