/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
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

#include <alps/ngs/observablewrappers.hpp>
// #ifdef ALPS_NGS_USE_NEW_ALEA
//     #include <alps/ngs/alea.hpp>
// #endif

namespace alps {

    namespace ngs {

        namespace detail {

            std::string ObservableWapper::getName() const {
                return _name;
            }

            uint32_t ObservableWapper::getBinnum() const {
                return _binnum;
            }

            std::string SignedObservableWapper::getSign() const {
                return _sign;
            }

        }
        
        //TODO
        alps::mcobservables & operator<< (alps::mcobservables & set, RealObservable const & obs) {
            set.create_RealObservable(obs.getName(), obs.getBinnum());
            return set;
        }

        // #ifdef ALPS_NGS_USE_NEW_ALEA
        //     alps::accumulator::accumulator_set & operator<< (alps::accumulator::accumulator_set & set, RealObservable const & obs) {
        //         using namespace alps::accumulator::tag;
                
        //         typedef accumulator::accumulator<double, accumulator::features<mean, error, max_num_binning> > accum_type;
        //         typedef accumulator::detail::accumulator_wrapper wrapper_type;
                
        //         set.insert(obs.getName(), boost::shared_ptr<wrapper_type>(new wrapper_type(accum_type(accumulator::bin_num = obs.getBinnum()))));

        //         return set;
        //     }
            
        //     alps::accumulator::accumulator_set & operator<< (alps::accumulator::accumulator_set & set, RealVectorObservable const & obs) {
        //         using namespace alps::accumulator::tag;
                
        //         typedef accumulator::accumulator<std::vector<double>, accumulator::features<mean, error, max_num_binning> > accum_type;
        //         typedef accumulator::detail::accumulator_wrapper wrapper_type;
                
        //         set.insert(obs.getName(), boost::shared_ptr<wrapper_type>(new wrapper_type(accum_type(accumulator::bin_num = obs.getBinnum()))));

        //         return set;
        //     }
            
        //     alps::accumulator::accumulator_set & operator<< (alps::accumulator::accumulator_set & set, SimpleRealObservable const & obs) {
        //         using namespace alps::accumulator::tag;
                
        //         typedef accumulator::accumulator<double, accumulator::features<mean, error> > accum_type;
        //         typedef accumulator::detail::accumulator_wrapper wrapper_type;
                
        //         set.insert(obs.getName(), boost::shared_ptr<wrapper_type>(new wrapper_type(accum_type(accumulator::bin_num = obs.getBinnum()))));
                
        //         return set;
        //     }
            
        //     alps::accumulator::accumulator_set & operator<< (alps::accumulator::accumulator_set & set, SimpleRealVectorObservable const & obs) {
        //         using namespace alps::accumulator::tag;
                
        //         typedef accumulator::accumulator<std::vector<double>, accumulator::features<mean, error> > accum_type;
        //         typedef accumulator::detail::accumulator_wrapper wrapper_type;
                
        //         set.insert(obs.getName(), boost::shared_ptr<wrapper_type>(new wrapper_type(accum_type(accumulator::bin_num = obs.getBinnum()))));
                
        //         return set;
        //     }
        // #endif        

        alps::mcobservables & operator<< (alps::mcobservables & set, RealVectorObservable const & obs) {
            set.create_RealVectorObservable(obs.getName(), obs.getBinnum());
            return set;
        }

        alps::mcobservables & operator<< (alps::mcobservables & set, SimpleRealObservable const & obs) {
            set.create_SimpleRealObservable(obs.getName());
            return set;
        }

        alps::mcobservables & operator<< (alps::mcobservables & set, SimpleRealVectorObservable const & obs) {
            set.create_SimpleRealVectorObservable(obs.getName());
            return set;
        }

        alps::mcobservables & operator<< (alps::mcobservables & set, SignedRealObservable const & obs) {
            set.create_SignedRealObservable(obs.getName(), obs.getSign(), obs.getBinnum());
            return set;
        }

        alps::mcobservables & operator<< (alps::mcobservables & set, SignedRealVectorObservable const & obs) {
            set.create_SignedRealVectorObservable(obs.getName(), obs.getSign(), obs.getBinnum());
            return set;
        }

        alps::mcobservables & operator<< (alps::mcobservables & set, SignedSimpleRealObservable const & obs) {
            set.create_SignedSimpleRealObservable(obs.getName(), obs.getSign());
            return set;
        }

        alps::mcobservables & operator<< (alps::mcobservables & set, SignedSimpleRealVectorObservable const & obs) {
            set.create_SignedSimpleRealVectorObservable(obs.getName(), obs.getSign());
            return set;
        }

        alps::mcobservables & operator<< (alps::mcobservables & set, RealTimeSeriesObservable const & obs) {
            set.create_RealTimeSeriesObservable(obs.getName());
            return set;
        }

        alps::mcobservables & operator<< (alps::mcobservables & set, RealVectorTimeSeriesObservable const & obs) {
          set.create_RealTimeSeriesObservable(obs.getName());
          return set;
        }
    };

}
