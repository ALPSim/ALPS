/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef LOOPER_EVALUATOR_H
#define LOOPER_EVALUATOR_H

#include <alps/scheduler.h>
#ifdef HAVE_PARAPACK_13
# include <alps/parapack/serial.h>
#else
# include <alps/parapack/worker.h>
#endif

namespace looper {

class abstract_evaluator : public alps::parapack::simple_evaluator {
public:
  virtual ~abstract_evaluator() {}
  virtual void pre_evaluate(alps::ObservableSet& m, alps::Parameters const&,
    alps::ObservableSet const& m_in) const = 0;
  virtual void evaluate(alps::scheduler::MCSimulation&, alps::Parameters const&,
    boost::filesystem::path const&) const = 0;
  virtual void evaluate(alps::ObservableSet& m, alps::Parameters const&,
    alps::ObservableSet const& m_in) const = 0;
  virtual void evaluate(alps::ObservableSet&) const {};
};

} // end namespace looper

#endif // LOOPER_EVALUATOR_H
