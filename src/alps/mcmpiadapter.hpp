/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
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

#include <alps/ngs/config.hpp>

#if !defined(ALPS_MCMPIADAPTER_HPP) && defined(ALPS_HAVE_MPI)
#define ALPS_MCMPIADAPTER_HPP

#include <alps/ngs.hpp>
#include <alps/ngs/boost_mpi.hpp>
#include <alps/check_schedule.hpp>

namespace alps {

    template<typename Base, typename ScheduleChecker = alps::check_schedule> class mcmpiadapter : public Base {

        public:

            mcmpiadapter(
                  typename Base::parameters_type const & parameters
                , boost::mpi::communicator const & comm
                , ScheduleChecker const & check = ScheduleChecker()
            )
                : Base(parameters, comm.rank())
                , communicator(comm)
                , schedule_checker(check)
                , clone(comm.rank())
#ifndef ALPS_NGS_USE_NEW_ALEA
                , binnumber(parameters["BINNUMBER"] | 128)
#endif
            {}

            double fraction_completed() const {
                return fraction;
            }

            bool run(boost::function<bool ()> const & stop_callback) {
                bool done = false, stopped = false;
                do {
                    this->update();
                    this->measure();
                    if (stopped || schedule_checker.pending()) {
                        stopped = stop_callback(); 
                        double local_fraction = stopped ? 1. : Base::fraction_completed();
                        schedule_checker.update(fraction = boost::mpi::all_reduce(communicator, local_fraction, std::plus<double>()));
                        done = fraction >= 1.;
                    }
                } while(!done);
                return !stopped;
            }

            typename Base::results_type collect_results() const {
                return collect_results(this->result_names());
            }

            typename Base::results_type collect_results(typename Base::result_names_type const & names) const {
                typename Base::results_type partial_results;
                for(typename Base::result_names_type::const_iterator it = names.begin(); it != names.end(); ++it) {
                    #ifdef ALPS_NGS_USE_NEW_ALEA
                        if (communicator.rank() == 0) {
                            if (this->measurements[*it].count()) {
                                typename Base::observable_collection_type::value_type merged = this->measurements[*it];
                                merged.collective_merge(communicator, 0);
                                partial_results.insert(*it, merged.result());
                            } else
                                partial_results.insert(*it, this->measurements[*it].result());
                        } else if (this->measurements[*it].count())
                            this->measurements[*it].collective_merge(communicator, 0);
                    #else
                        alps::mcresult result(this->measurements[*it]); // TODO: use Base::collect_results
                        if (result.count())
                            partial_results.insert(*it, result.reduce(communicator, binnumber));
                        else
                            partial_results.insert(*it, result);
                    #endif
                }
                return partial_results;
            }

        protected:

            boost::mpi::communicator communicator;

            ScheduleChecker schedule_checker;
            double fraction;
            int clone;
#ifndef ALPS_NGS_USE_NEW_ALEA
            std::size_t binnumber;
#endif
    };
}

#endif
