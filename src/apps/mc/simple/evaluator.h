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

  void evaluate(alps::alea::result_set& obs) const {
    if (obs.has("Inverse Temperature") && obs.has("Number of Sites") &&
        obs.has("Energy") && obs.has("Energy^2")) {
      const alps::alea::observable_result beta = obs["Inverse Temperature"];
      const alps::alea::observable_result n = obs["Number of Sites"];
      const alps::alea::observable_result ene = obs["Energy"];
      const alps::alea::observable_result ene2 = obs["Energy^2"];
      if (beta.count() && n.count() && ene.count() && ene2.count())
        obs.insert("Specific Heat", beta.mean() * beta.mean() * (ene2 - ene * ene) / n.mean());
    }
    if (obs.has("Magnetization Density^2") && obs.has("Magnetization Density^4")) {
      const alps::alea::observable_result m2 = obs["Magnetization Density^2"];
      const alps::alea::observable_result m4 = obs["Magnetization Density^4"];
      if (m2.count() && m4.count())
        obs.insert("Binder Ratio of Magnetization", m2 * m2 / m4);
    }
    if (obs.has("Magnetization Density Z^2") && obs.has("Magnetization Density Z^4")) {
      const alps::alea::observable_result m2 = obs["Magnetization Density Z^2"];
      const alps::alea::observable_result m4 = obs["Magnetization Density Z^4"];
      if (m2.count() && m4.count())
        obs.insert("Binder Ratio of Magnetization Z", m2 * m2 / m4);
    }
    if (obs.has("Magnetization Density X^2") && obs.has("Magnetization Density X^4")) {
      const alps::alea::observable_result m2 = obs["Magnetization Density X^2"];
      const alps::alea::observable_result m4 = obs["Magnetization Density X^4"];
      if (m2.count() && m4.count())
        obs.insert("Binder Ratio of Magnetization X", m2 * m2 / m4);
    }
    if (obs.has("Inverse Temperature")) obs.erase("Inverse Temperature");
    if (obs.has("Temperature")) obs.erase("Temperature");
  }
};

#endif // ALPS_SIMPLE_MC_EVALUATOR_H
