/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef LOOPER_EVALUATOR_IMPL_H
#define LOOPER_EVALUATOR_IMPL_H

#include "evaluator.h"
#include "measurement.h"
#include <boost/foreach.hpp>

namespace looper {

template<typename MEASUREMENT_SET>
class evaluator : public abstract_evaluator {
public:
  typedef typename measurement<MEASUREMENT_SET>::type measurement_t;
  evaluator() {}
  evaluator(alps::Parameters const& p) : params(p) {}
  void pre_evaluate(alps::ObservableSet& m, alps::Parameters const& p,
    alps::ObservableSet const& m_in) const {
    pre_evaluator_selector<measurement_t>::pre_evaluate(m, p, m_in);
  }
  void evaluate(alps::scheduler::MCSimulation& sim, alps::Parameters const& p,
    boost::filesystem::path const&) const {
    alps::ObservableSet m;
    energy_evaluator::evaluate(m, sim.get_measurements());
    evaluator_selector<measurement_t>::evaluate(m, p, sim.get_measurements());
    BOOST_FOREACH(alps::ObservableSet::iterator::value_type const& v, m)
      sim.addObservable(*(v.second));
  }
  void evaluate(alps::ObservableSet& m, alps::Parameters const& p,
    alps::ObservableSet const& m_in) const {
    energy_evaluator::evaluate(m, m_in);
    evaluator_selector<measurement_t>::evaluate(m, p, m_in);
  }
  void evaluate(alps::ObservableSet& m) const {
    energy_evaluator::evaluate(m, m);
    evaluator_selector<measurement_t>::evaluate(m, params, m);
  }
private:
  alps::Parameters params;
};

} // end namespace looper

#endif // LOOPER_EVALUATOR_IMPL_H
