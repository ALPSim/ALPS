/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS/simple-mc
*
* Copyright (C) 1997-2015 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef ALPS_SIMPLE_MC_EVALUATOR_H
#define ALPS_SIMPLE_MC_EVALUATOR_H

#include <alps/parapack/worker.h>

class evaluator : public alps::parapack::simple_evaluator {
public:
  evaluator(alps::Parameters const&) {}
  virtual ~evaluator() {}

  void evaluate(alps::ObservableSet& obs) const {
    if (obs.has("Inverse Temperature") && obs.has("Number of Sites") &&
        obs.has("Energy") && obs.has("Energy^2")) {
      alps::RealObsevaluator beta = obs["Inverse Temperature"];
      alps::RealObsevaluator n = obs["Number of Sites"];
      alps::RealObsevaluator ene = obs["Energy"];
      alps::RealObsevaluator ene2 = obs["Energy^2"];
      if (beta.count() && n.count() && ene.count() && ene2.count()) {
        alps::RealObsevaluator c("Specific Heat");
        c = beta.mean() * beta.mean() * (ene2 - ene * ene) / n.mean();
        obs.addObservable(c);
      }
    }
    if (obs.has("Magnetization Density^2") && obs.has("Magnetization Density^4")) {
      alps::RealObsevaluator m2 = obs["Magnetization Density^2"];
      alps::RealObsevaluator m4 = obs["Magnetization Density^4"];
      if (m2.count() && m4.count()) {
        alps::RealObsevaluator binder("Binder Ratio of Magnetization");
        binder = m2 * m2 / m4;
        obs.addObservable(binder);
      }
    }
    if (obs.has("Magnetization Density Z^2") && obs.has("Magnetization Density Z^4")) {
      alps::RealObsevaluator m2 = obs["Magnetization Density Z^2"];
      alps::RealObsevaluator m4 = obs["Magnetization Density Z^4"];
      if (m2.count() && m4.count()) {
        alps::RealObsevaluator binder("Binder Ratio of Magnetization Z");
        binder = m2 * m2 / m4;
        obs.addObservable(binder);
      }
    }
    if (obs.has("Magnetization Density X^2") && obs.has("Magnetization Density X^4")) {
      alps::RealObsevaluator m2 = obs["Magnetization Density X^2"];
      alps::RealObsevaluator m4 = obs["Magnetization Density X^4"];
      if (m2.count() && m4.count()) {
        alps::RealObsevaluator binder("Binder Ratio of Magnetization X");
        binder = m2 * m2 / m4;
        obs.addObservable(binder);
      }
    }
    if (obs.has("Inverse Temperature")) obs.erase("Inverse Temperature");
    if (obs.has("Temperature")) obs.erase("Temperature");
  }
};

#endif // ALPS_SIMPLE_MC_EVALUATOR_H
