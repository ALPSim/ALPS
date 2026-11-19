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

namespace looper {

template<typename MEASUREMENT_SET>
class evaluator : public abstract_evaluator {
public:
  typedef typename measurement<MEASUREMENT_SET>::type measurement_t;
  evaluator() {}
  evaluator(alps::Parameters const& p) : params(p) {}
  void pre_evaluate(alps::alea::result_set& m, alps::Parameters const& p,
    alps::alea::result_set const& m_in) const {
    pre_evaluator_selector<measurement_t>::pre_evaluate(m, p, m_in);
  }
  void evaluate(alps::scheduler::MCSimulation& sim, alps::Parameters const& p,
    std::filesystem::path const&) const {
    alps::alea::result_set m;
    const alps::alea::result_set m_in = sim.get_measurements();
    energy_evaluator::evaluate(m, m_in);
    evaluator_selector<measurement_t>::evaluate(m, p, m_in);
    for (std::string const& name : m.names())
      sim.add_result(name, m[name]);
  }
  void evaluate(alps::alea::result_set& m, alps::Parameters const& p,
    alps::alea::result_set const& m_in) const {
    energy_evaluator::evaluate(m, m_in);
    evaluator_selector<measurement_t>::evaluate(m, p, m_in);
  }
  void evaluate(alps::alea::result_set& m) const {
    const alps::alea::result_set m_in(m);
    energy_evaluator::evaluate(m, m_in);
    evaluator_selector<measurement_t>::evaluate(m, params, m_in);
  }
private:
  alps::Parameters params;
};

} // end namespace looper

#endif // LOOPER_EVALUATOR_IMPL_H
