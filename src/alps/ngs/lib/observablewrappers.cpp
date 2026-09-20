/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/observablewrappers.hpp>

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
