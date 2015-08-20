/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS/simple-mc
*
* Copyright (C) 1997-2015 by Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
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
    if (obs.has("Magnetization^2") && obs.has("Magnetization^4")) {
      alps::RealObsevaluator m2 = obs["Magnetization^2"];
      alps::RealObsevaluator m4 = obs["Magnetization^4"];
      if (m2.count() && m4.count()) {
        alps::RealObsevaluator binder("Binder Ratio of Magnetization");
        binder = m2 * m2 / m4;
        obs.addObservable(binder);
      }
    }
    if (obs.has("Magnetization Z^2") && obs.has("Magnetization Z^4")) {
      alps::RealObsevaluator m2 = obs["Magnetization Z^2"];
      alps::RealObsevaluator m4 = obs["Magnetization Z^4"];
      if (m2.count() && m4.count()) {
        alps::RealObsevaluator binder("Binder Ratio of Magnetization Z");
        binder = m2 * m2 / m4;
        obs.addObservable(binder);
      }
    }
    if (obs.has("Magnetization X^2") && obs.has("Magnetization X^4")) {
      alps::RealObsevaluator m2 = obs["Magnetization X^2"];
      alps::RealObsevaluator m4 = obs["Magnetization X^4"];
      if (m2.count() && m4.count()) {
        alps::RealObsevaluator binder("Binder Ratio of Magnetization X");
        binder = m2 * m2 / m4;
        obs.addObservable(binder);
      }
    }
  }
};

#endif // ALPS_SIMPLE_MC_EVALUATOR_H
